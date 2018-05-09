#In this script, I'll use MethylKit to analyze differences in methylation between C. virginica gonad samples.

#### SET WORKING DIRECTORY ####
getwd()
setwd("../2018-04-27-Bismark/") #Set working directory as bismark folder

#### INSTALL PACKAGES ####
install.packages("devtools") #Install the devtools package
library(devtools) #Load devtools

source("https://bioconductor.org/biocLite.R") #Source package from bioconductor
biocLite("methylKit") #Install methylkit
library(methylKit) #Load methylkit

install_github("al2na/methylKit", build_vignettes = FALSE, repos=BiocInstaller::biocinstallRepos(), dependencies = TRUE) #Install more methylKit options

#### PROCESS METHYLATION DATA ####

analysisFiles <- list ("zr2096_1_s1_R1_bismark_bt2_pe.deduplicated.sorted.bam",
                       "zr2096_2_s1_R1_bismark_bt2_pe.deduplicated.sorted.bam",
                       "zr2096_3_s1_R1_bismark_bt2_pe.deduplicated.sorted.bam",
                       "zr2096_4_s1_R1_bismark_bt2_pe.deduplicated.sorted.bam",
                       "zr2096_5_s1_R1_bismark_bt2_pe.deduplicated.sorted.bam",
                       "zr2096_6_s1_R1_bismark_bt2_pe.deduplicated.sorted.bam",
                       "zr2096_7_s1_R1_bismark_bt2_pe.deduplicated.sorted.bam",
                       "zr2096_8_s1_R1_bismark_bt2_pe.deduplicated.sorted.bam",
                       "zr2096_9_s1_R1_bismark_bt2_pe.deduplicated.sorted.bam",
                       "zr2096_10_s1_R1_bismark_bt2_pe.deduplicated.sorted.bam") #Put all .bam files into a list for analysis
sample.IDs <- list("1", "2", "3", "4", "5", "6", "7", "8", "9", "10") #Create list of sample IDs
treatmentSpecification <- c(rep(0, times = 5), rep(1, times = 5))
processedFiles <- processBismarkAln(location = analysisFiles, sample.id = sample.IDs, assembly = "v3", read.context = "CpG", mincov = 1, treatment = treatmentSpecification) #Process files for CpG meetehylation. First 5 files were from ambient conditionds, and the second from high pCO2 conditions.

#### ANALYZE METHYLATION DATA ####

####NEED TO CHANGE WORKING DIRECTORY THEN RUN FILES. MOVE GRAPHS FROM BISMARK FOLDER ####

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
  jpeg(filename = fileName$actualFileName[i], height = 1000, width = 1000) #Save file with designated name
  getCoverageStats(processedFiles[[i]], plot = TRUE, both.strands = FALSE) #Get CpG coverage information
  dev.off() #Turn off plotting device
} #Plot and save CpG coverage information

methylationInformation <- unite(processedFiles)

?unite

getCorrelation(methylationInformation,plot=TRUE)

clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)

getCorrelation(meth,plot=TRUE)


clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)



hc = clusterSamples(meth, dist="correlation", method="ward", plot=FALSE)

PCASamples(meth, screeplot=TRUE)

PCASamples(meth)

myDiff=calculateDiffMeth(meth)

myDiff25p=getMethylDiff(myDiff,difference=25,qvalue=0.01)

myDiff50p <- getMethylDiff(myDiff,difference=50,qvalue=0.01)

write.table(myDiff50p, file = "analyses/myDiff50p.tab")

View(myDiff50p)