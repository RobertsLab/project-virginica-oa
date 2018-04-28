# `analyses` Subdirectory Structure

This subdirectory hosts the output of all analyses regarding *C. virginica* and ocean acidification. Every folder houses all necessary files for a discrete analysis step in a larger pipeline. Folders follow this naming convention:

**Year-Month-Date-AnalysisStep** *(ex. 2018-01-23-MBDSeq-Labwork)*

The following pipelines are represented in this subdirectory:

**[2018-01-23-MBDSeq-Labwork](https://github.com/RobertsLab/project-virginica-oa/blob/master/analyses/2018-01-23-MBDSeq-Labwork/)**:

- [Labwork Calculations](https://github.com/RobertsLab/project-virginica-oa/blob/master/analyses/2018-01-23-MBDSeq-Labwork/2018-01-23-Virginica-MBDSeq-Labwork-Calculations.xlsx)

**Gonad Methylation Analyses**:

- [2018-04-26-Gonad-Methylation-FastQC](https://github.com/RobertsLab/project-virginica-oa/tree/master/analyses/2018-04-26-Gonad-Methylation-FastQC): Quality-checks for sequences
  - [MultiQC Report](https://github.com/RobertsLab/project-virginica-oa/blob/master/analyses/2018-04-26-Gonad-Methylation-FastQC/multiqc_report.html)
  - [MultiQC Data](https://github.com/RobertsLab/project-virginica-oa/tree/master/analyses/2018-04-26-Gonad-Methylation-FastQC/multiqc_data)
- [2018-04-27-Gonad-Methylation-Bismark](https://github.com/RobertsLab/project-virginica-oa/tree/master/analyses/2018-04-27-Bismark): Sequence alignment output from `bismark`. Genome preparation outputs are included in `.gitignore`.
