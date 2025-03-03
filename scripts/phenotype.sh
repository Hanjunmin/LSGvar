# snv
x=$1 
y=$2 
ref=$3
bcftools merge -m none --force-samples "${x}/sortsnv.vcf.gz" "${y}/sortsnv.vcf.gz" > mergesnv.vcf
file="mergesnv.vcf" 

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
file="mergesnvend.vcf" 
if [ -f "$file.gz" ]; then
  rm "$file.gz"
fi

bgzip $file
bcftools sort $file".gz" -o "sort"$file".gz"
bcftools index -t "sort"$file".gz"

bcftools merge -m none --force-samples "${x}/sortindel.vcf.gz" "${y}/sortindel.vcf.gz" |bgzip  >mergeindel.vcf.gz
file="mergeindel.vcf" 
bcftools sort $file".gz" -o "sort"$file".gz"
bcftools index -t "sort"$file".gz"
truvari collapse -s 0 -S 50 -p 0.7 -P 0.7 -r 10 -i "sort"$file".gz" -o "${file}_merge.vcf" -c "${file}_collapsed.vcf" -f  $ref
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
file="mergeindelend.vcf" 
if [ -f "$file.gz" ]; then
  rm "$file.gz"
fi

bgzip $file
bcftools sort $file".gz" -o "sort"$file".gz"
bcftools index -t "sort"$file".gz"

# SV  
bcftools merge -m none --force-samples "${x}/sortSV.vcf.gz" "${y}/sortSV.vcf.gz" |bgzip >mergeSV.vcf.gz
file="mergeSV.vcf" 
bcftools sort $file".gz" -o "sort"$file".gz"
bcftools index -t "sort"$file".gz"
truvari collapse -i "sort"$file".gz" -o mergesv.vcf -c collapsedsv.vcf -f $ref
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
file="mergesv_end.vcf" 
if [ -f "$file.gz" ]; then
  rm "$file.gz"
fi
bgzip $file
bcftools sort $file".gz" -o "sort"$file".gz"
bcftools index -t "sort"$file".gz"

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

gunzip -f sortLSGvar.vcf.gz
awk '{print $0 "\t" NR}' sortLSGvar.vcf >sortLSGvar.num.vcf
bedtools intersect -a sortLSGvar.num.vcf -b h1_h2intersec.bed |uniq > overlapped.vcf
less -S overlapped.vcf  |awk '{print $11}' |uniq >overlap.num
end=$(less sortLSGvar.vcf |wc -l)
seq 1 ${end} >all.num
sort all.num overlap.num| uniq -u |sort -n >sediff.num
paste <(less -S  overlapped.vcf |awk 'OFS="\t"{print $1,$2,$3,$4,$5,$6,$7,$8,$9}') <(less -S  overlapped.vcf |cut -f 10 |sed 's/\./0/g') >over.vcf
awk 'NR==FNR {lines[$1]; next} FNR in lines' sediff.num sortLSGvar.vcf >sed.vcf
cat sed.vcf over.vcf >LSGvarall.vcf
file="LSGvarall.vcf" 
if [ -f "$file.gz" ]; then
  rm "$file.gz"
fi
bgzip $file
bcftools sort $file".gz" -o "sort"$file".gz"
bcftools index -t "sort"$file".gz"
cp  "sort"$file".gz" ../
