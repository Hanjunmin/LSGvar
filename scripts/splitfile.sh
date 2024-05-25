folder=$1
mkdir ${folder}
#sed -i '1,5s/^#[^#]/#&/' $2
# 38:36
# 19:40
awk -F'\t' '{
    if($0 ~ /^#/){  
        print $0
    } else{
            OFS="\t";
            $11=(length($4)-length($5) >0 ? length($4)-length($5) : length($5)-length($4))
            print $1, $2, $3, $4, $5, $6, $7, $8, $9,$10, $11
        }
}' $2   > "${folder}/count.vcf"  #"${nowdic}/macigarsdr.vcf"

less -S "${folder}/count.vcf" |awk -F'\t' '{
   if($0 ~ /^#/){ 
        print $0
    }else{
            OFS="\t";
            if($11>=50){
                 print $1, $2, $3, $4, $5, $6, $7, $8, $9,$10
                }else{
                    next
                }
        }
}' >"${folder}/SV.vcf"

less -S "${folder}/count.vcf"  |awk -F'\t' '{
    if($0 ~ /^#/){ 
        print $0
    }else{
            OFS="\t";
            if($11<50 && $11>0){
                 print $1, $2, $3, $4, $5, $6, $7, $8, $9,$10
                }else{
                    next
                }
        }
}' >"${folder}/indel.vcf"

less -S "${folder}/count.vcf"  |awk -F'\t' '{
     if($0 ~ /^#/){ 
        print $0
    }else{
            OFS="\t";
            if($11==0){
                 print $1, $2, $3, $4, $5, $6, $7, $8, $9,$10
                }else{
                    next
                }
        }
}' >"${folder}/snv.vcf"

## SV
cd ${folder}
file="SV.vcf" #hg002benend.vcf ourend.vcf giab.vcf ourpaend.vcf output_file.vcf
rm $file".gz"
bgzip $file
bcftools sort $file".gz" -o "sort"$file".gz"
bcftools index -t "sort"$file".gz"

file="snv.vcf" #hg002benend.vcf ourend.vcf giab.vcf ourpaend.vcf output_file.vcf
rm $file".gz"
bgzip $file
bcftools sort $file".gz" -o "sort"$file".gz"
bcftools index -t "sort"$file".gz"

file="indel.vcf" #hg002benend.vcf ourend.vcf giab.vcf ourpaend.vcf output_file.vcf
rm $file".gz"
bgzip $file
bcftools sort $file".gz" -o "sort"$file".gz"
bcftools index -t "sort"$file".gz"
