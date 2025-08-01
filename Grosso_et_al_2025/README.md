# Causal relationships between gut microbiome and hundreds of age-related traits: evidence of a replicable effect on ApoM protein levels

**Release date: 03 February 2025**

**Date of last update: 11 June 2025**

Scripts used for Mendelian randomization analysis in the paper *Grosso et al. (2025)* (**running on R version 4.4.1**).

The pipeline is divided in 3 folders containing the codes used to perform the main Mendelian randomization analysis. 
1. The files in the "***GWAS_QC***" folder were first used to correct missing values using the 1000 Genomes European reference panel and other corrections related to the GWASs used for the project specifically. Then there are files used to perform clumping.
2. The "***MR_analysis***" folder contains the codes for the main MR analysis, sensitivity analyses, plots and FDR correction.
3. In the folder "***Power_analysis_and_figures***" there is the Rmd file to perform the power analysis of significant results and the respective replication analyses.

<< **Use of AI statement:**
We acknowledge the use of artificial intelligence tools, namely ChatGPT 4.0, to improve code efficiency, particularly to reduce the development time of clear and accurate result tables. >>


***QUESTIONS?*** 
For questions regarding the code please contact Federica Grosso (federicagrosso@cnr.it) or Serena Sanna (serena.sanna@cnr.it)
