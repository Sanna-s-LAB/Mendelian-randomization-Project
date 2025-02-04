# Causal relationships between gut microbiome and age-related traits
### Date released: 03/02/2025
Scripts used for Mendelian randomization analysis in the paper *Grosso et al. (2025)* (**running on R version 4.4.1**).

**$MedR\chi v$ doi: https://doi.org/10.1101/2025.02.03.25321568**

The pipeline is divided in 3 folders containing the codes used to perform the main Mendelian randomization analysis. 
1. The files in the "***GWAS_QC***" folder were first used to correct missing values using the 1000 Genomes European reference panel and other corrections related to the GWASs used for the project specifically. Then there are files used to perform clumping.
2. The "***MR_analysis***" folder contains the codes for the main MR analysis, sensitivity analyses, plots and FDR correction.
3. In the folder "***Power_analysis***" there is the Rmd file to perform the power analysis of significant results and the respective replication analyses.

<< **Use of AI statement:**
We acknowledge the use of artificial intelligence tools, namely ChatGPT 4.0, to improve code efficiency, particularly to reduce the development time of clear and accurate result tables. >>


***QUESTIONS?*** 
For questions regarding the code please contact Federica Grosso (federicagrosso@cnr.it) or Serena Sanna (serena.sanna@cnr.it)
