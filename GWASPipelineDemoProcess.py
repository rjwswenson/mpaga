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

#get data
hl.utils.get_1kg('pdata/')

#Create MatrixTable DF with downsample 1000 Genomes Data
mt = hl.read_matrix_table('pdata/1kg.mt')

#Prepare Table 
table = (hl.import_table('pdata/1kg_annotations.txt', impute=True)
         .key_by('Sample'))
         
#Add Annotations 
mt = mt.annotate_cols(pheno = table[mt.s])
#Aggregate
mt.aggregate_cols(hl.agg.counter(mt.pheno.SuperPopulation))

#Count SNPs
snp_counts = mt.aggregate_rows(hl.agg.counter(hl.Struct(ref=mt.alleles[0], alt=mt.alleles[1])))

from collections import Counter
counts = Counter(snp_counts)
counts.most_common()

#Prepare Sample QC 
mt = hl.sample_qc(mt)

#Filtering QC 
mt = mt.filter_cols((mt.sample_qc.dp_stats.mean >= 4) & (mt.sample_qc.call_rate >= 0.97))


#Genotypic QC
call_rate = mt.aggregate_entries(hl.agg.fraction(hl.is_defined(mt.GT)))

ab = mt.AD[1] / hl.sum(mt.AD)
filter_condition_ab = ((mt.GT.is_hom_ref() & (ab <= 0.1)) |
                        (mt.GT.is_het() & (ab >= 0.25) & (ab <= 0.75)) |
                        (mt.GT.is_hom_var() & (ab >= 0.9)))

mt = mt.filter_entries(filter_condition_ab)

#Post QC Call Rate
post_qc_call_rate = mt.aggregate_entries(hl.agg.fraction(hl.is_defined(mt.GT)))


#Variant QC
mt = hl.variant_qc(mt)


#Perform GWAS
#Cutoff of 1%
mt = mt.filter_rows(mt.variant_qc.AF[1] > 0.01)
#Not so far from Hardy-Weinberg equilibrium
mt = mt.filter_rows(mt.variant_qc.p_value_hwe > 1e-6)
#Identify number of downsamples and variants

#GWAS with Phenotype of Interest: Caffeine Consumption
gwas = hl.linear_regression_rows(y=mt.pheno.CaffeineConsumption,
                                 x=mt.GT.n_alt_alleles(),
                                 covariates=[1.0])



#Confounded!
eigenvalues, pcs, _ = hl.hwe_normalized_pca(mt.GT)


#Annotate
mt = mt.annotate_cols(scores = pcs[mt.s].scores)


#Rerun GWAS with linear regression, controlling for sample sex 
gwas = hl.linear_regression_rows(
    y=mt.pheno.CaffeineConsumption,
    x=mt.GT.n_alt_alleles(),
    covariates=[1.0, mt.pheno.isFemale, mt.scores[0], mt.scores[1], mt.scores[2]])




#Rare Variants Analysis
entries = mt.entries()
results = (entries.group_by(pop = entries.pheno.SuperPopulation, chromosome = entries.locus.contig)
      .aggregate(n_het = hl.agg.count_where(entries.GT.is_het())))


entries = entries.annotate(maf_bin = hl.cond(entries.info.AF[0]<0.01, "< 1%",
                             hl.cond(entries.info.AF[0]<0.05, "1%-5%", ">5%")))

results2 = (entries.group_by(af_bin = entries.maf_bin, purple_hair = entries.pheno.PurpleHair)
      .aggregate(mean_gq = hl.agg.stats(entries.GQ).mean,
                 mean_dp = hl.agg.stats(entries.DP).mean))


#Store our MatrixTable mt onto HDFS
mt.write('pout/gwas-completed.mt')
#Convert our Hail MatrixTable into a VCF "VCF 4.2 Specification" File
hl.export_vcf(mt, 'pout/gwas-completed.vcf.bgz') 















