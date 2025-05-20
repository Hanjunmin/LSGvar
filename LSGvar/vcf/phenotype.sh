#!/bin/bash

hap1=$1  ## hap1 folder
hap2=$2  ## hap2 folder
ref=$3   ## reference genome
variants=${4:-all}

cd temp || exit 1

process_snv=false
process_indel=false
process_inv=false

if [ "$variants" = "all" ]; then
  process_snv=true
  process_indel=true
  process_inv=true
else
  IFS=',' read -r -a variant_arr <<< "$variants"
  for v in "${variant_arr[@]}"; do
    case "$v" in
      snv) process_snv=true ;;
      ins) process_indel=true ;;
      del) process_indel=true ;;
      inv) process_inv=true ;;
      *)
        echo "Warning: unknown variant type '$v' ignored."
        ;;
    esac
  done
fi

# --------- merge snv -----------
if $process_snv; then
  bcftools merge -m none --force-samples ../"${hap1}/sortsnv.vcf.gz" ../"${hap2}/sortsnv.vcf.gz" > mergesnv.vcf

  awk -F '\t' '{
      if ($0 ~ /^##/) {
          print $0;
      } else if ($0 ~ /^#CHROM/) {
          print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample";
      } else {
          OFS="\t";
          if ($10 == "1|0") {
              $12 = ($11 == "1|0") ? "1|1" : "1|.";
          } else {
              $12 = ".|1";
          }
          $10 = "";
          $11 = "";
          print $1, $2, $3, $4, $5, $6, $7, $8, $9, $12;
      }
  }' mergesnv.vcf > mergesnvend.vcf

  [ -f mergesnvend.vcf.gz ] && rm mergesnvend.vcf.gz
  bgzip mergesnvend.vcf
  bcftools sort mergesnvend.vcf.gz -o sortmergesnvend.vcf.gz
  bcftools index -t sortmergesnvend.vcf.gz
fi

# --------- merge indel -----------
if $process_indel; then
  bcftools merge -m none --force-samples ../"${hap1}/sortindel.vcf.gz" ../"${hap2}/sortindel.vcf.gz" | bgzip > mergeindel.vcf.gz
  bcftools sort mergeindel.vcf.gz -o sortmergeindel.vcf.gz
  bcftools index -t sortmergeindel.vcf.gz

  truvari collapse -s 0 -S 50 -p 0.7 -P 0.7 -r 10 -i sortmergeindel.vcf.gz -o mergeindel_merge.vcf -c mergeindel_collapsed.vcf -f $ref

  awk -F '\t' '{
      if ($0 ~ /^##/) {
          print $0;
      } else if ($0 ~ /^#CHROM/) {
          print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample";
      } else {
          OFS="\t";
          if ($10 == "1|0") {
              $12 = ($11 == "1|0") ? "1|1" : "1|.";
          } else {
              $12 = ".|1";
          }
          $10 = "";
          $11 = "";
          print $1, $2, $3, $4, $5, $6, $7, $8, $9, $12;
      }
  }' mergeindel_merge.vcf > mergeindelend.vcf

  [ -f mergeindelend.vcf.gz ] && rm mergeindelend.vcf.gz
  bgzip mergeindelend.vcf
  bcftools sort mergeindelend.vcf.gz -o sortmergeindelend.vcf.gz
  bcftools index -t sortmergeindelend.vcf.gz
fi

# --------- merge inv -----------
if $process_inv; then
  bcftools merge -m none --force-samples ../"${hap1}/sortSV.vcf.gz" ../"${hap2}/sortSV.vcf.gz" | bgzip > mergeSV.vcf.gz
  bcftools sort mergeSV.vcf.gz -o sortmergeSV.vcf.gz
  bcftools index -t sortmergeSV.vcf.gz

  bcftools norm -m -both sortmergeSV.vcf.gz | bgzip > mergeSV_split.vcf.gz
  bcftools sort mergeSV_split.vcf.gz -o sortmergeSV_split.vcf.gz
  tabix -p vcf sortmergeSV_split.vcf.gz

  truvari collapse -i sortmergeSV_split.vcf.gz -o mergesv.vcf -c collapsedsv.vcf -f $ref

  awk -F '\t' '{
      if ($0 ~ /^##/) {
          print $0;
      } else if ($0 ~ /^#CHROM/) {
          print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample";
      } else {
          OFS="\t";
          if ($10 == "1|0") {
              $12 = ($11 == "1|0") ? "1|1" : "1|.";
          } else {
              $12 = ".|1";
          }
          $10 = "";
          $11 = "";
          print $1, $2, $3, $4, $5, $6, $7, $8, $9, $12;
      }
  }' mergesv.vcf > mergesv_end.vcf

  [ -f mergesv_end.vcf.gz ] && rm mergesv_end.vcf.gz
  bgzip mergesv_end.vcf
  bcftools sort mergesv_end.vcf.gz -o sortmergesv_end.vcf.gz
  bcftools index -t sortmergesv_end.vcf.gz
fi

# --------- merge all variants -----------

cat_header() {
  zcat "$1" | grep '^##'
  echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample"
}

cat_vcfs() {
  local files=("$@")
  local header_printed=false

  for f in "${files[@]}"; do
    if [ -f "$f" ]; then
      if ! $header_printed; then
        cat_header "$f"
        header_printed=true
      fi
      zcat "$f" | grep -v '^##' | grep -v '^#CHROM'
    fi
  done
}

vcf_files=()
$process_inv && vcf_files+=("sortmergesv_end.vcf.gz")
$process_indel && vcf_files+=("sortmergeindelend.vcf.gz")
$process_snv && vcf_files+=("sortmergesnvend.vcf.gz")

cat_vcfs "${vcf_files[@]}" > LSGvar.vcf

[ -f LSGvar.vcf.gz ] && rm LSGvar.vcf.gz
bgzip LSGvar.vcf
bcftools sort LSGvar.vcf.gz -o sortLSGvar.vcf.gz
bcftools index -t sortLSGvar.vcf.gz

gunzip -f sortLSGvar.vcf.gz
awk '{print $0 "\t" NR}' sortLSGvar.vcf > sortLSGvar.num.vcf
bedtools intersect -a sortLSGvar.num.vcf -b ../results/integrate/h1_h2intersec.bed | uniq > overlapped.vcf
awk '{print $11}' overlapped.vcf | uniq > overlap.num
end=$(wc -l < sortLSGvar.vcf)
seq 1 $end > all.num
sort all.num overlap.num | uniq -u | sort -n > sediff.num
paste <(awk 'OFS="\t"{print $1,$2,$3,$4,$5,$6,$7,$8,$9}' overlapped.vcf) <(cut -f 10 overlapped.vcf | sed 's/\./0/g') > over.vcf
awk 'NR==FNR {lines[$1]; next} FNR in lines' sediff.num sortLSGvar.vcf > sediff.vcf
cat sediff.vcf over.vcf > LSGvarall.vcf

[ -f LSGvarall.vcf.gz ] && rm LSGvarall.vcf.gz
bgzip LSGvarall.vcf
bcftools sort LSGvarall.vcf.gz -o sortLSGvarall.vcf.gz
bcftools index -t sortLSGvarall.vcf.gz

cp sortLSGvarall.vcf.gz ../results/

cd ..
rm -r temp/
