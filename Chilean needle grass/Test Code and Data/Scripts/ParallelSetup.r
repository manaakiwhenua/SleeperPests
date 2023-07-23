rm(list = ls())

######################################################
###Important: 
###1)make sure using latest version of R or 
###devtools might not install
###2) install all packages used in foreach loop to the default
###R library location. "Workers" in the cluster do not recognise custom paths 
######################################################

######################################
###Install required packages for 
###for running CNG implementation of INApest
######################################
install.packages("doParallel")
install.packages("abind")
install.packages("devtools")
install.packages("igraph")
install.packages("tidyverse")
install.packages("actuar")
devtools::install_github("GarrettLab/INA")
