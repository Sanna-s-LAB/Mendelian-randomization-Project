lipids=(TG HDC LDC nonHDC TC)

#sex=FEMALE
#statins=1

sex=$1
statins=$2


cd /groups/umcg-lifelines/tmp01/projects/ov23_0790/results_regenie/${sex}_${statins}/ 
for l in ${lipids[@]}
do
    echo $l
    rm all_chr_${sex}_${statins}_step2_${l}.regenie.txt
    rm all_chr_${sex}_${statins}_step2_${l}.regenie.txt.gz
    
    for f in *_${sex}_${statins}_step2_${l}.regenie_srt.txt.gz
    do
        zcat $f | awk -v lip=${l} 'BEGIN {FS=OFS="\t"}; {print lip, $0}'  >> all_chr_${sex}_${statins}_step2_${l}.regenie.txt
    done
    sort -k15,15g all_chr_${sex}_${statins}_step2_${l}.regenie.txt | gzip -c > all_chr_${sex}_${statins}_step2_${l}.regenie.txt.gz
    rm all_chr_${sex}_${statins}_step2_${l}.regenie.txt
done

