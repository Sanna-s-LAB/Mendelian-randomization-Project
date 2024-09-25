1. Taking clumping results and entire GWAS of outcome, start with ***run_analysis.sh*** to run analysis of **MR.R** (with a "for" loop each of 37 GWAS of exposure vs 1 outcome in this case) and ***FDR_correction.R*** to correct IVW p-values 
2. This analysis includes:
   - Harmonization
   - Main MR analysis with IVW, MR-Egger, Weighted median
   - Pleiotropy
   - Heterogeneity
   - Leave-one-out with leave-one-out plot
   - Scatter plot
   - MR-PRESSO (using for plots the functions modified from package *TwoSampleMR* and putted in the "***MR_PRESSO_plots.R***" script)
   all this results are saved in files in the same folder, named with the accession number of exposure and outcome
   
