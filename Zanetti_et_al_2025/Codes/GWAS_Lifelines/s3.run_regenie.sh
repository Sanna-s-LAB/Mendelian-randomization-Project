#!/bin/bash
#SBATCH --job-name=regenie
#SBATCH --output=logs/run_GWAS.out
#SBATCH --error=logs/run_GWAS.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=30gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L


ml regenie

snplist=/groups/umcg-lifelines/tmp01/projects/ov23_0790/data/genotypes/gsa_typed.txt

chr=$1
sex=$2
statins=$3
d=/groups/umcg-lifelines/tmp01/projects/ov23_0790/

echo "chr=$chr, sex=$sex, statins=$statins"

#sex=FEMALE
#statins=0
pheno_file=${d}/data/phenotypes/pheno_adj_${sex}_${statins}.txt
geno_file=${d}/data/genotypes/merged/chr${chr}.filtered
covar_file=${d}/data/phenotypes/genotype_batch_covariate.txt

res_file=${d}/results_regenie/${sex}_${statins}/${chr}_${sex}_${statins}
mkdir -p ${d}/results_regenie/${sex}_${statins}/

# Step 1
regenie \
    --step 1 \
    --bed ${geno_file} \
    --extract ${snplist} \
    --phenoFile ${pheno_file} \
    --covarFile ${covar_file} \
    --bsize 200 \
    --out ${res_file}_step1 \
    --loocv


# Step 2
PRED=${res_file}_step1_pred.list
regenie \
    --step 2 \
    --bed  ${geno_file} \
    --phenoFile ${pheno_file} \
    --covarFile ${covar_file} \
    --pred ${PRED} \
    --bsize 200 \
    --out ${res_file}_step2

rm ${res_file}_step1*


## Postprocessing
ml PythonPlus/3.10.4-GCCcore-11.3.0-v23.01.1

for f in ${d}/results_regenie/${sex}_${statins}/${chr}_${sex}_${statins}*.regenie
do
  awk 'BEGIN {OFS="\t"}; {if (NR == 1) {$13="PVAL"; print} else {$13=10^-$12; gsub(/:[ACGT]_[ACGT]/, "", $3); print $0}}' $f | \
    python ~/scripts/umcg_scripts/misc/add_rs_by_position.py stdin /groups/umcg-llnext/tmp01/umcg-dzhernakova/resources/dbsnp_137.b37.vcf.gz 2 | \
    sort -k14,14g | gzip -c > ${f}_srt.txt.gz
  rm $f
done