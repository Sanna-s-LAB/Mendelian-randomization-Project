#
# 1. Filter SNPs by imputation quality (for each batch)
#

ml PLINK

plink2 \
  --vcf /groups/umcg-lifelines/tmp01/releases/gsa_imputed/v2/vcf_files/${CHR}.GSAr2.vcf.gz \
  --extract-if-info "INFO > 0.4" \
  --make-bed \
  --out /groups/umcg-lifelines/tmp01/projects/ov23_0790/data/genotypes/UGLI1/chr${CHR}.r2_0.4
  
plink2 \
  --vcf /groups/umcg-lifelines/tmp01/projects/ov23_0790/data/genotypes/UGLI2/tmp_vcfs/${CHR}.UGLI2_r2.vcf.gz \
  --extract-if-info "INFO > 0.4" \
  --make-bed \
  --out /groups/umcg-lifelines/tmp01/projects/ov23_0790/data/genotypes/UGLI2//chr${CHR}.r2_0.4


#
# 2. Change SNP ids to chr:pos (for each batch)
#
cd /groups/umcg-lifelines/tmp01/projects/ov23_0790/data/genotypes/UGLI1/
for f in *bim
do
    echo $f
    mv $f ${f}.orig
    awk 'BEGIN {FS=OFS="\t"}; {$2 = $1 ":" $4 ":" $5 "_" $6; print }' ${f}.orig > ${f}
done

cd /groups/umcg-lifelines/tmp01/projects/ov23_0790/data/genotypes/UGLI2/
for f in *bim
do
    echo $f
    mv $f ${f}.orig
    awk 'BEGIN {FS=OFS="\t"}; {$2 = $1 ":" $4 ":" $5 "_" $6; print }' ${f}.orig > ${f}
done


#
# 3. Merge 2 genotype batches and filter
#
cd /groups/umcg-lifelines/tmp01/projects/ov23_0790/data/genotypes/merged
ml PLINK/1.9-beta6-20190617

for chr in `seq 1 22`
do
    plink --bfile ../UGLI1/chr${chr}.r2_0.4 --bmerge ../UGLI2/chr${chr}.r2_0.4 --make-bed --out chr${chr}

    plink --bfile chr${chr} --maf 0.05 --geno 0.05 --hwe 1e-6 --make-bed --out chr${chr}.filtered

done

#
# 4. Get genotyped (non-imputed) SNPs for regenie
#

rm gsa_typed.txt
for f in /groups/umcg-lifelines/prm03/releases/gsa_genotypes/v2/VCF/*vcf.gz
do 
    echo $f
    zcat $f | awk 'BEGIN {FS=OFS="\t"}; {if (!/^#/) print $1 ":" $2 ":" $5 "_" $4}' >> gsa_typed.txt
done

rm affymetrix_typed.txt
for f in UGLI2/tmp_vcfs/*
do 
  echo $f
  bcftools view -G -H -i 'INFO/TYPED==1'  $f | cut -f1-5 >> affymetrix_typed.txt
done