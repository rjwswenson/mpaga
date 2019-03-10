#@Ryan J.W. Swenson
#rswenson@cloudera.com 
#CLDR Hail 1KG Genome GWAS Demo

#Cloudera Hail 1KG Genome GWAS Demo
from pyspark import SparkContext, SparkConf
from pyspark.sql import SQLContext
import pyspark.sql.functions as sf
from pyspark.sql.functions import udf
conf = SparkConf().setAppName("app")
#Pauses for Context Load
sc = SparkContext(conf=conf)
sqlContext = SQLContext(sc)


#Import and Initialize Hail 
import hail as hl
print ('Initializing Hail...')
hl.init(sc)

#Load Bokeh IO Visualization
from bokeh.io import output_notebook, show
from pprint import pprint
output_notebook()

#Load Downsample 1000 Genomes Data From Google Cloud Storage S3 Bucket
#hl.utils.get_1kg('1000genomesdata/')

mt = hl.read_matrix_table('1000genomesdata/1kg.mt')


mt.rows().select().show(5)
mt.row_key.show(5)
mt.s.show(5)

#Prepare Table 
table = (hl.import_table('1000genomesdata/1kg_annotations.txt', impute=True)
         .key_by('Sample'))
         
#Show Table
table.show(width=100)

#Add Annotations 
mt = mt.annotate_cols(pheno = table[mt.s])
#Aggregate
mt.aggregate_cols(hl.agg.counter(mt.pheno.SuperPopulation))

#Count SNPs
snp_counts = mt.aggregate_rows(hl.agg.counter(hl.Struct(ref=mt.alleles[0], alt=mt.alleles[1])))
pprint(snp_counts)


from collections import Counter
counts = Counter(snp_counts)
counts.most_common()

#Plot DP Histogram
p = hl.plot.histogram(mt.DP, range=(0,30), bins=30, title='DP Histogram', legend='DP')
show(p)

#Describe Columns of MatrixTable
mt.col.describe()

#Prepare Sample QC 
mt = hl.sample_qc(mt)

#Schema should now have been expanded to include QC
mt.col.describe()

#Sample QC Call Rates (.88 to 1) Plot
p = hl.plot.histogram(mt.sample_qc.call_rate, range=(.88,1), legend='Call Rate')
show(p)

#Mean Sample QC (Phreds 10 - 70) Plot
p = hl.plot.histogram(mt.sample_qc.gq_stats.mean, range=(10,70), legend='Mean Sample GQ')
show(p)

#Show Correlation if exist
p = hl.plot.scatter(mt.sample_qc.dp_stats.mean, mt.sample_qc.call_rate, xlabel='Mean DP', ylabel='Call Rate')
show(p)

#Filtering QC 
mt = mt.filter_cols((mt.sample_qc.dp_stats.mean >= 4) & (mt.sample_qc.call_rate >= 0.97))
print('After filter, %d/284 samples remain.' % mt.count_cols())

#Genotypic QC
call_rate = mt.aggregate_entries(hl.agg.fraction(hl.is_defined(mt.GT)))
print('before genotype QC, call rate is %.3f' % call_rate)


# if we find a genotype called homozygous reference with >10% alternate reads, 
#a genotype called homozygous alternate with >10% reference reads, 
#or a genotype called heterozygote without a ref / alt balance near 1:1,
#it is likely to be an error
ab = mt.AD[1] / hl.sum(mt.AD)
filter_condition_ab = ((mt.GT.is_hom_ref() & (ab <= 0.1)) |
                        (mt.GT.is_het() & (ab >= 0.25) & (ab <= 0.75)) |
                        (mt.GT.is_hom_var() & (ab >= 0.9)))

mt = mt.filter_entries(filter_condition_ab)

#Post QC Call Rate
post_qc_call_rate = mt.aggregate_entries(hl.agg.fraction(hl.is_defined(mt.GT)))
print('post QC call rate is %.3f' % post_qc_call_rate)

#Variant QC
mt = hl.variant_qc(mt)
mt.row.describe()

#Perform GWAS
#Cutoff of 1%
mt = mt.filter_rows(mt.variant_qc.AF[1] > 0.01)
#Not so far from Hardy-Weinberg equilibrium
mt = mt.filter_rows(mt.variant_qc.p_value_hwe > 1e-6)
#Identify number of downsamples and variants
print('Samples: %d  Variants: %d' % (mt.count_cols(), mt.count_rows()))


#GWAS with Phenotype of Interest: Caffeine Consumption
gwas = hl.linear_regression_rows(y=mt.pheno.CaffeineConsumption,
                                 x=mt.GT.n_alt_alleles(),
                                 covariates=[1.0])
gwas.row.describe()
#Perform Manhatten Plot
p = hl.plot.manhattan(gwas.p_value)
show(p)
#Perform Quantile-Quantile Plot
p = hl.plot.qq(gwas.p_value)
show(p)

#Confounded!
eigenvalues, pcs, _ = hl.hwe_normalized_pca(mt.GT)
pprint(eigenvalues)
pcs.show(5, width=100)

#Annotate
mt = mt.annotate_cols(scores = pcs[mt.s].scores)

#Plot
p = hl.plot.scatter(mt.scores[0],
                    mt.scores[1],
                    label=mt.pheno.SuperPopulation,
                    title='PCA', xlabel='PC1', ylabel='PC2')
show(p)

#Rerun GWAS with linear regression, controlling for sample sex 
gwas = hl.linear_regression_rows(
    y=mt.pheno.CaffeineConsumption,
    x=mt.GT.n_alt_alleles(),
    covariates=[1.0, mt.pheno.isFemale, mt.scores[0], mt.scores[1], mt.scores[2]])

#Updated Quantile-Quantile Plot
p = hl.plot.qq(gwas.p_value)
show(p)

#Updated Manhatten Plot
p = hl.plot.manhattan(gwas.p_value)
show(p)

#Rare Variants Analysis
entries = mt.entries()
results = (entries.group_by(pop = entries.pheno.SuperPopulation, chromosome = entries.locus.contig)
      .aggregate(n_het = hl.agg.count_where(entries.GT.is_het())))


entries = entries.annotate(maf_bin = hl.cond(entries.info.AF[0]<0.01, "< 1%",
                             hl.cond(entries.info.AF[0]<0.05, "1%-5%", ">5%")))

results2 = (entries.group_by(af_bin = entries.maf_bin, purple_hair = entries.pheno.PurpleHair)
      .aggregate(mean_gq = hl.agg.stats(entries.GQ).mean,
                 mean_dp = hl.agg.stats(entries.DP).mean))

results2.show()

#Store our MatrixTable mt onto HDFS
mt.write('hailv2_outputs/gwas.mt')

#Convert our Hail MatrixTable into a VCF "VCF 4.2 Specification" File
hl.export_vcf(mt, 'hailv2_outputs/gwas.vcf.bgz') 















