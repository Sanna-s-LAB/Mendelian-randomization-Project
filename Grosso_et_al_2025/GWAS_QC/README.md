These files are used to to the standard quality control of GWAS files pre-MR and clumping.
- First delete the rows where there were $-4<\beta<4$ (the first few rows) and correct a wrong base pair position. (Using ***rm_beta.sh*** to run ***RM_beta.R***) $\to$ Create a folder **GWAS_new**, where to put new GWASs created
- Then since there were NAs in the microbiome files use the reference file of the 1000Genomes EUR from which you took the columns "*chr:position*" and rsid. Merge with all GWAS using chr:pos key. In this way it can be obteined a new rsid column. When the original column was NA, it replaced it with the new rsid. Then the code creates a new column where 
    - 0 = rsid are equal
    - 1 = one of them is NA
    - 2 = rsids are different

This is done running ***ref_col.sh*** file to run ***merge.R***.
    
- Finally, it creates a "rsid_to_use" column with all the rsids of the microbiome when they are there and are the same as those of the reference, those of the 1KGP reference when they are different or when the microbiome has NA → ***final_column.sh*** which put final files in the **GWAS_final** folder

- Use the function ***run_clumping.sh*** to run ***Clumping.R*** to do clumping (in this case with $p=5\cdot 10^{-6}$, but usually with $p=5\cdot 10^{-8}$)
