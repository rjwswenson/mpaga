#!/bin/bash
add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/' &>>install.log
apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 &>>install.log
apt update &>>install.log
apt install r-base-dev -y &>> install.log
apt install r-base -y &>> install.log
