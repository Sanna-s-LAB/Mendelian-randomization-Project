These are the scripts used for the Mendelian Randomization (MR) analysis and for presenting the results:
1. First, we run the **run_MR_analysis.sh** Bash script, which executes the **MR_analysis.R** code. This performs the MR analysis for all pairs of proteins and lipid outcomes, separately for women and men, and generates the corresponding results tables for each pair. The same script can be adapted for different types of analyses: cis+trans MR, cis-only MR, trans-only MR, bidirectional MR, or Lifelines, by changing the column names accordingly.
2. Next, the **Tables.R** script is used to produce the final overall result tables with all results together. This script relies on several functions defined in **process_tables_functions.R**.
3. Finally, the **Forest_plots.R** script generates forest plots to visually represent the final results for both women and men together.

