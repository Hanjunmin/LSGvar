
# snv
less -S $4 |awk 'OFS="\t"{print $1,$3,$4}'|bedtools sort |bedtools merge -i - >h1_paf.bed
less -S $5 |awk 'OFS="\t"{print $1,$3,$4}'|bedtools sort |bedtools merge -i - >h2_paf.bed
bedtools intersect -a h1_paf.bed -b h2_paf.bed  |bedtools sort >h1_h2intersec.bed


x=$1 #x="/home/jmhan/SDR/HG002/GRCH37/ma"  x="/home/jmhan/SDR/HG002/ma"
y=$2 #y="/home/jmhan/SDR/HG002/GRCH37/pa"  y="/home/jmhan/SDR/HG002/pa"  
#ref="/home/jmhan/SDR/genome/GRCh38/GCF_000001405.26_GRCh38_genomic.chrNames.fna"
#ref="/home/jmhan/SDR/HG002/GRCH37/hg19.fasta"
ref=$3
bcftools merge -m none --force-samples "${x}/sortsnv.vcf.gz" "${y}/sortsnv.vcf.gz" > mergesnv.vcf
file="mergesnv.vcf" #hg002benend.vcf ourend.vcf giab.vcf ourpaend.vcf output_file.vcf
# hg38:37
# 19:42
awk -F '\t' '{
    if ($0 ~ /^##/) {
        print $0; 
    } else if ($0 ~ /^#CHROM/) {
        print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002";
    } else {
        OFS="\t";  
        if ($10 == "1|0") {
            $12 = ($11 == "1|0") ? "1|1" : "1|.";
        } else {
            $12 = ".|1";
        }
        $10 = "";
        $11 = "";
        print $1, $2, $3, $4, $5, $6, $7, $8, $9, $12
    }
}' mergesnv.vcf > mergesnvend.vcf
file="mergesnvend.vcf" #hg002benend.vcf ourend.vcf giab.vcf ourpaend.vcf output_file.vcf
if [ -f "$file.gz" ]; then
  rm "$file.gz"
fi

bgzip $file
bcftools sort $file".gz" -o "sort"$file".gz"
bcftools index -t "sort"$file".gz"



bcftools merge -m none --force-samples "${x}/sortindel.vcf.gz" "${y}/sortindel.vcf.gz" |bgzip  >mergeindel.vcf.gz
file="mergeindel.vcf" #hg002benend.vcf ourend.vcf giab.vcf ourpaend.vcf output_file.vcf
bcftools sort $file".gz" -o "sort"$file".gz"
bcftools index -t "sort"$file".gz"
/share/home/zhanglab/user/yangchentao/miniconda3/bin/truvari collapse -s 0 -S 50 -p 0.7 -P 0.7 -r 10 -i "sort"$file".gz" -o "${file}_merge.vcf" -c "${file}_collapsed.vcf" -f  $ref
awk -F '\t' '{
    if ($0 ~ /^##/) {
        print $0;  
    } else if ($0 ~ /^#CHROM/) {
        print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002";
    } else {
        OFS="\t";  
        if ($10 == "1|0") {
            $12 = ($11 == "1|0") ? "1|1" : "1|.";
        } else {
            $12 = ".|1";
        }
        $10 = "";
        $11 = "";
        print $1, $2, $3, $4, $5, $6, $7, $8, $9, $12
    }
}' "${file}_merge.vcf" > mergeindelend.vcf
file="mergeindelend.vcf" #hg002benend.vcf ourend.vcf giab.vcf ourpaend.vcf output_file.vcf
if [ -f "$file.gz" ]; then
  rm "$file.gz"
fi

bgzip $file
bcftools sort $file".gz" -o "sort"$file".gz"
bcftools index -t "sort"$file".gz"
#nohup truvari bench -s 0 --sizemax 100000000 -b /home/jmhan/SDR/HG002/CIGAR/HG002process/sorthg002benindel.vcf.gz -c "sort"$file".gz" --reference /home/jmhan/SDR/genome/GRCh38/GCF_000001405.26_GRCh38_genomic.chrNames.fna --includebed /home/jmhan/SDR/HG002/genome/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed -o delme &





# SV  ## hg19
bcftools merge -m none --force-samples "${x}/sortSV.vcf.gz" "${y}/sortSV.vcf.gz" |bgzip >mergeSV.vcf.gz
file="mergeSV.vcf" #hg002benend.vcf ourend.vcf giab.vcf ourpaend.vcf output_file.vcf
bcftools sort $file".gz" -o "sort"$file".gz"
bcftools index -t "sort"$file".gz"
#conda activate truvari4
/share/home/zhanglab/user/yangchentao/miniconda3/bin/truvari collapse -i "sort"$file".gz" -o mergesv.vcf -c collapsedsv.vcf -f $ref
awk -F '\t' '{
    if ($0 ~ /^##/) {
        print $0; 
    } else if ($0 ~ /^#CHROM/) {
        print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002";
    } else {
        OFS="\t";  
        if ($10 == "1|0") {
            $12 = ($11 == "1|0") ? "1|1" : "1|.";
        } else {
            $12 = ".|1";
        }
        $10 = "";
        $11 = "";
        print $1, $2, $3, $4, $5, $6, $7, $8, $9,$12
    }
}' mergesv.vcf > mergesv_end.vcf
file="mergesv_end.vcf" #hg002benend.vcf ourend.vcf giab.vcf ourpaend.vcf output_file.vcf
rm $file".gz"
bgzip $file
bcftools sort $file".gz" -o "sort"$file".gz"
bcftools index -t "sort"$file".gz"
#nohup truvari bench -b /home/jmhan/SDR/HG002/genome/hg002sv_12745/sortbenv0.6.vcf.gz -c "sort"$file".gz" --reference /home/jmhan/SDR/HG002/GRCH37/hg19.fasta --includebed /home/jmhan/SDR/HG002/genome/v0.6consistent.bed -o delme &

## integrate

cat <(zcat sortmergesv_end.vcf.gz|awk -F'\t' '{if($0 ~ /^#/){print $0}else{OFS="\t"
        print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}}') <(zcat sortmergeindelend.vcf.gz |awk -F'\t' '{if($0 ~ /^#/){next}else{OFS="\t"
        print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}}') <(zcat sortmergesnvend.vcf.gz |awk -F'\t' '{if($0 ~ /^#/){next}else{OFS="\t"
        print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}}') > LSGvar.vcf
file="LSGvar.vcf" #hg002benend.vcf ourend.vcf giab.vcf ourpaend.vcf output_file.vcf
if [ -f "$file.gz" ]; then
  rm "$file.gz"
fi

bgzip $file
bcftools sort $file".gz" -o "sort"$file".gz"
bcftools index -t "sort"$file".gz"


gunzip sortLSGvar.vcf.gz
awk '{print $0 "\t" NR}' sortLSGvar.vcf >sortLSGvar.num.vcf
bedtools intersect -a sortLSGvar.num.vcf -b h1_h2intersec.bed |uniq > overlapped.vcf
less -S overlapped.vcf  |awk '{print $11}' |uniq >overlap.num
end=$(less sortLSGvar.vcf |wc -l)
seq 1 ${end} >all.num
sort all.num overlap.num| uniq -u |sort -n >sediff.num
paste <(less -S  overlapped.vcf |awk 'OFS="\t"{print $1,$2,$3,$4,$5,$6,$7,$8,$9}') <(less -S  overlapped.vcf |cut -f 10 |sed 's/\./0/g') >over.vcf
awk 'NR==FNR {lines[$1]; next} FNR in lines' sediff.num sortLSGvar.vcf >sed.vcf
cat sed.vcf over.vcf >LSGvarall.vcf
file="LSGvarall.vcf" #hg002benend.vcf ourend.vcf giab.vcf ourpaend.vcf output_file.vcf
if [ -f "$file.gz" ]; then
  rm "$file.gz"
fi
bgzip $file
bcftools sort $file".gz" -o "sort"$file".gz"
bcftools index -t "sort"$file".gz"
cp "sort"$file".gz" ../


## grch38
# cat <(zcat /home/jmhan/SDR/HG002/integrate/sortmergesnvend.vcf.gz) <(zcat /home/jmhan/SDR/HG002/integrate/sortmergeindelend.vcf.gz | awk -F'\t' '{
#     if ($0 ~ /^#/) {
#         next
#     } else {
#         OFS="\t"
#         print $1, $2, $3, $4, $5, $6, $7, ".", $9, $10
#     }
# }'
# ) > Lasv.vcf
# file="Lasv.vcf" #hg002benend.vcf ourend.vcf giab.vcf ourpaend.vcf output_file.vcf
# rm $file".gz"
# bgzip $file
# bcftools sort $file".gz" -o "sort"$file".gz"
# bcftools index -t "sort"$file".gz"
#nohup truvari bench -b /home/jmhan/SDR/HG002/genome/hg002sv_12745/sortbenv0.6.vcf.gz -c "sort"$file".gz" --reference /home/jmhan/SDR/HG002/GRCH37/hg19.fasta --includebed /home/jmhan/SDR/HG002/genome/v0.6consistent.bed -o end &
#singularity shell /home/jmhan/singularity/hap.py_latest.sif
#/opt/hap.py/bin/som.py sorthg002benend.vcf.gz sorthg002benend.vcf.gz -o de -r /home/jmhan/SDR/genome/GRCh38/GCF_000001405.26_GRCh38_genomic.chrNames.fna
#singularity exec /home/jmhan/SDR/HG002/integrate/hap.py/hap.py/hap.py_latest.sif /opt/hap.py/bin/hap.py   hg002_38_splitben.vcf       hg002_38_splitben.vcf      -f /home/jmhan/SDR/HG002/genome/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed -r  /home/jmhan/SDR/genome/GRCh38/GCF_000001405.26_GRCh38_genomic.chrNames.fna   -o test/test


# less -S SV.vcf |awk -F'\t' '{
#     if(NR<36){
#         print $0
#     }else if (NR==36){
#         print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002";
#         }else{
#             OFS="\t";
           
#                 print $1, $2, ".", $4, $5, ".", ".", ".", $9,$10
#         }
# }' >paSV.vcf




# filein="/home/jmhan/SDR/HG002/GRCH37/sortmergesv_end.vcf.gz"
# bedtools intersect <(less -S  ${filein} |awk -F'\t' '{
#     if(NR < 46){
#         next
#     } else{
# OFS="\t";
#             $11=(length($4)-length($5) >0 ? length($4)-length($5) : length($5)-length($4))
#             print $1, $2, $2+$11
#     }
# }') <(less -S fn.vcf.gz |awk 'NR>=111{print $0}' |less -S


# )
