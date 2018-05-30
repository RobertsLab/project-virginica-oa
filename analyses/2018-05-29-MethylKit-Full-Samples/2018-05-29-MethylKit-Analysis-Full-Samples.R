#In this script, I'll use MethylKit to analyze differences in methylation between C. virginica gonad samples.

#### SET WORKING DIRECTORY ####
#getwd()
#setwd("../2018-05-22-Bismark-Full-Samples//") #Set working directory as bismark folder

#### INSTALL PACKAGES ####
install.packages("devtools") #Install the devtools package
library(devtools) #Load devtools

source("https://bioconductor.org/biocLite.R") #Source package from bioconductor
biocLite("methylKit") #Install methylkit
install_github("al2na/methylKit", build_vignettes = FALSE, repos=BiocInstaller::biocinstallRepos(), dependencies = TRUE) #Install more methylKit options
library(methylKit) #Load methylkit

#### PROCESS METHYLATION DATA ####

analysisFiles <- list ("zr2096_1_dedup.sorted.bam",
                       "zr2096_2_dedup.sorted.bam",
                       "zr2096_3_dedup.sorted.bam",
                       "zr2096_4_dedup.sorted.bam",
                       "zr2096_5_dedup.sorted.bam",
                       "zr2096_6_dedup.sorted.bam",
                       "zr2096_7_dedup.sorted.bam",
                       "zr2096_8_dedup.sorted.bam",
                       "zr2096_9_dedup.sorted.bam",
                       "zr2096_10_dedup.sorted.bam") #Put all .bam files into a list for analysis
sample.IDs <- list("1", "2", "3", "4", "5", "6", "7", "8", "9", "10") #Create list of sample IDs
treatmentSpecification <- c(rep(0, times = 5), rep(1, times = 5))
processedFiles <- processBismarkAln(location = analysisFiles, sample.id = sample.IDs, assembly = "v3", read.context = "CpG", mincov = 1, treatment = treatmentSpecification) #Process files for CpG meetehylation. First 5 files were from ambient conditions, and the second from high pCO2 conditions.

#### ANALYZE METHYLATION DATA ####

setwd("../2018-05-29-MethylKit-Full-Samples/") #Change working direcotry to the MethylKit folder

nFiles <- length(sample.IDs) #Count number of samples
fileName <- data.frame("nameBase" = rep("2018-05-30-Percent-CpG-Methylation", times = nFiles),
                       "nameBase2" = rep("2018-05-30-Percent-CpG-Coverage", times = nFiles),
                       "sample.ID" = 1:10) #Create new dataframe for filenames
head(fileName) #Confirm dataframe creation
fileName$actualFileName <- paste(fileName$nameBase, "-Sample", fileName$sample.ID, ".jpeg", sep = "") #Create a new column for the full filename
fileName$actualFileName2 <- paste(fileName$nameBase2, "-Sample", fileName$sample.ID, ".jpeg", sep = "") #Create a new column for the full filename
head(fileName) #Confirm column creation

for(i in 1:nFiles) { #For each data file
  jpeg(filename = fileName$actualFileName[i], height = 1000, width = 1000) #Save file with designated name
  getMethylationStats(processedFiles[[i]], plot = TRUE, both.strands = FALSE) #Get %CpG methylation information
  dev.off() #Turn off plotting device
} #Plot and save %CpG methylation information

for(i in 1:nFiles) { #For each data file
  jpeg(filename = fileName$actualFileName2[i], height = 1000, width = 1000) #Save file with designated name
  getCoverageStats(processedFiles[[i]], plot = TRUE, both.strands = FALSE) #Get CpG coverage information
  dev.off() #Turn off plotting device
} #Plot and save CpG coverage information

methylationInformation <- unite(processedFiles) #Combine all processed files into a single table

#jpeg(filename = "2018-05-30-Full-Sample-Pearson-Correlation-Plot.jpeg", height = 1000, width = 1000) #Save file with designated name
getCorrelation(methylationInformation, plot = TRUE) #Understand correlation between methylation patterns in different samples
#dev.off()

#jpeg(filename = "2018-05-30-Full-Sample-CpG-Methylation-Clustering.jpeg", height = 1000, width = 1000) #Save file with designated name
clusterSamples(methylationInformation, dist = "correlation", method = "ward", plot = TRUE) #Cluster samples based on correlation coefficients
#dev.off()

clusteringInformation <- clusterSamples(methylationInformation, dist = "correlation", method = "ward", plot = FALSE) #Save cluster information as a new object

#jpeg(filename = "2018-05-30-Full-Sample-Methylation-PCA.jpeg", height = 1000, width = 1000) #Save file with designated name
PCASamples(methylationInformation) #Run a PCA analysis on percent methylation for all samples
#dev.off() #Turn off plotting device

#jpeg(filename = "2018-05-30-Full-Sample-Methylation-Screeplot.jpeg", height = 1000, width = 1000) #Save file with designated name
PCASamples(methylationInformation, screeplot = TRUE) #Run the PCA analysis and plot variances against PC number in a screeplot
#dev.off()

differentialMethylationStats <- calculateDiffMeth(methylationInformation) #Calculate differential methylation statistics based on treatment indication from processBismarkAln
diffMethStats25 <- getMethylDiff(differentialMethylationStats, difference = 25, qvalue = 0.01) #Identify loci that are at least 25% different. Q-value is the FDR used for p-value corrections.
diffMethStats50 <- getMethylDiff(differentialMethylationStats, difference = 50, qvalue = 0.01) #Identify loci that are at least 50% different
head(diffMethStats50) #Confirm creation
#write.csv(diffMethStats50, "2018-05-30-Differentially-Methylated-Loci-50.csv") #Save table as .csv

#### SAVE DMLs AS A BED FILE ####

library(readr) #Load package
library(tidyverse) #Load package

diffMethStats50 <- read.csv("../2018-05-29-MethylKit-Full-Samples/2018-05-30-Differentially-Methylated-Loci-50.csv") #Import data
head(diffMethStats50) #Confirm import

DMLPlus05302018 <- filter(diffMethStats50, strand == "+") %>% mutate(start = start -1, end = end + 1) %>% select(chr, start, end, strand, meth.diff) #Save + strand of DMLs as a new object
DMLMinus05302018 <- filter(diffMethStats50, strand == "-") %>% mutate(start = start -2) %>% select(chr, start, end, strand, meth.diff) #Save - strand of DMLs as a new object
DML05302018 <- bind_rows(DMLPlus05302018, DMLMinus05302018) %>% arrange(chr, start) %>% mutate_if(is.numeric, as.integer) #Join + and - strand information to be saved as a BED file, and avoid writing information in scientific notation
write_delim(DML05302018, "2018-05-30-DML-Locations.bed",  delim = '\t', col_names = FALSE) #Save data as a BED file