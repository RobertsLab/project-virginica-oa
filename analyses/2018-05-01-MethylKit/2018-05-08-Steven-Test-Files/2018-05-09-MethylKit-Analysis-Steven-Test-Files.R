#In this script, I'll test the methylKit pipeline with 2 files Steven generated. This will allow me to figure out why I cannot unite my files. If I can unite Steven's files, then I know there is an issue with the files I generated. If I cannot unite Steven's files, then it could be a software issue.

#### SET WORKING DIRECTORY ####
getwd()
setwd("2018-05-08-Steven-Test-Files/") #Set working directory as test file folder

#### INSTALL PACKAGES ####
install.packages("devtools") #Install the devtools package
library(devtools) #Load devtools

source("https://bioconductor.org/biocLite.R") #Source package from bioconductor
biocLite("methylKit") #Install methylkit
install_github("al2na/methylKit", build_vignettes = FALSE, repos=BiocInstaller::biocinstallRepos(), dependencies = TRUE) #Install more methylKit options
library(methylKit) #Load methylkit

#### PROCESS METHYLATION DATA ####

analysisFiles <- list ("zr2096_1_dedup.sorted.bam",
                       "zr2096_10_dedup.sorted.bam") #Put all .bam files into a list for analysis
sample.IDs <- list("1", "10") #Create list of sample IDs
treatmentSpecification <- c(0, 1)
processedFiles <- processBismarkAln(location = analysisFiles, sample.id = sample.IDs, assembly = "v3", read.context = "CpG", mincov = 1, treatment = treatmentSpecification) #Process files for CpG meetehylation. First 5 files were from ambient conditionds, and the second from high pCO2 conditions.

#### ANALYZE METHYLATION DATA ####

nFiles <- length(sample.IDs) #Count number of samples
fileName <- data.frame("nameBase" = rep("2018-05-08-Percent-CpG-Methylation", times = nFiles),
                       "nameBase2" = rep("2018-05-08-Percent-CpG-Coverage", times = nFiles),
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

getCorrelation(methylationInformation, plot = TRUE) #Understand correlation between methylation patterns in different samples
clusterSamples(methylationInformation, dist = "correlation", method = "ward", plot = TRUE) #Cluster samples based on correlation coefficients
clusteringInformation <- clusterSamples(methylationInformation, dist = "correlation", method = "ward", plot = FALSE) #Save cluster information as a new object

PCASamples(methylationInformation) #Run a PCA analysis on percent methylation for all samples
PCASamples(methylationInformation, screeplot = TRUE) #Run the PCA analysis and plot variances against PC number in a screeplot

differentialMethylationStats <- calculateDiffMeth(methylationInformation) #Calculate differential methylation statistics based on treatment indication from processBismarkAln
diffMethStats25 <- getMethylDiff(differentialMethylationStats, difference = 25, qvalue = 0.01) #Identify loci that are at least 25% different
diffMethStats50 <- getMethylDiff(differentialMethylationStats,difference=50,qvalue=0.01) #Identify loci that are at least 50% different