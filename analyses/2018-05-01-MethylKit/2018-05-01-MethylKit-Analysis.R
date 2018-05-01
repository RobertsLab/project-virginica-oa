#In this script, I'll use MethylKit to analyze differences in methylation between C. virginica gonad samples.

#### SET WORKING DIRECTORY ####
getwd()
setwd("../") #Set working directory as analyses folder

#### INSTALL METHYLKIT ####
source("https://bioconductor.org/biocLite.R") #Source package from bioconductor
biocLite("methylKit")
install.packages("data.table")