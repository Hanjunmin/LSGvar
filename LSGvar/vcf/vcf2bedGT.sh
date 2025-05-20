#!/bin/bash

file=$1  ## sv file
hap1_sv_file=$2  ## hap1 sv bed
hap2_sv_file=$3  ## hap2 sv bed
output_bed=$4    ## output bed
variants=${5:-all}

# # Specific variants merge
process_snv=false
process_ins=false
process_del=false
process_inv=false
process_trans=false
process_sdr=false
process_dup=false
process_highdup=false
process_complex=false

if [ "$variants" = "all" ]; then
  process_snv=true
  process_ins=true
  process_del=true
  process_inv=true
  process_trans=true
  process_sdr=true
  process_dup=true
  process_highdup=true
  process_complex=true
else
  IFS=',' read -r -a variant_arr <<< "$variants"
  for v in "${variant_arr[@]}"; do
    case "$v" in
      snv) process_snv=true ;;
      ins) process_ins=true ;;
      del) process_del=true ;;
      inv) process_inv=true ;;
      trans) process_trans=true ;;
      sdr) process_sdr=true ;;
      dup) process_dup=true ;;
      highdup) process_highdup=true ;;
      complex) process_complex=true ;;
      *)
        echo "Warning: unknown variant type '$v' ignored."
        ;;
    esac
  done
fi

mkdir -p results/integrate

if $process_snv; then
  paste <(less -S $file | grep 'SNV' | awk 'OFS="\t"{print $1,$2,$3,$2+1,"SNV",1,".",$10}' | less -S) \
        <(less -S $file | grep 'SNV' | awk -F'TIG_REGION=' '{print $2}' | less -S | awk -F',' '{print $1}') \
        > results/integrate/SNV.bed
fi

if $process_ins; then
  paste <(less -S $file | grep 'INS' | awk 'OFS="\t"{print $1,$2,$3}') \
        <(less -S $file | grep 'INS' | awk -F'[-\t]' 'OFS="\t"{print $2+1,"INS",$6,"."}' | less -S) \
        <(less -S $file | grep 'INS' | cut -f10) \
        <(less -S $file | grep 'INS' | awk -F'TIG_REGION=' '{print $2}' | less -S | awk -F',' '{print $1}') \
        > results/integrate/INS.bed
fi

if $process_del; then
  paste <(less -S $file | grep 'DEL' | awk 'OFS="\t"{print $1,$2,$3}') \
        <(less -S $file | grep 'DEL' | awk -F'[-\t]' 'OFS="\t"{print $2+$6,"DEL",$6,"."}' | less -S) \
        <(less -S $file | grep 'DEL' | cut -f10) \
        <(less -S $file | grep 'DEL' | awk -F'TIG_REGION=' '{print $2}' | less -S | awk -F',' '{print $1}') \
        > results/integrate/DEL.bed
fi

if $process_inv; then
  paste <(less -S $file | grep 'INV' | awk 'OFS="\t"{print $1,$2,$3}') \
        <(less -S $file | grep 'INV' | awk -F'[-\t]' 'OFS="\t"{print $2+$6,"INV",$6,"."}' | less -S) \
        <(less -S $file | grep 'INV' | cut -f10) \
        <(less -S $file | grep 'INV' | awk -F'TIG_REGION=' '{print $2}' | less -S | awk -F',' '{print $1}') \
        > results/integrate/INV.bed
fi

# merge
bed_files=""
$process_snv && bed_files+="results/integrate/SNV.bed "
$process_ins && bed_files+="results/integrate/INS.bed "
$process_del && bed_files+="results/integrate/DEL.bed "
$process_inv && bed_files+="results/integrate/INV.bed "
echo "SV split done!"
cat $bed_files > results/integrate/mid.bed

awk 'OFS="\t"{if ($8 == "0|1" || $8 == ".|1") print $1,$2,$3,$4,$5,$6,"hap2",$8,$9; else if ($8 == "1|0" || $8 == "1|.") print $1,$2,$3,$4,$5,$6,"hap1",$8,$9; else if ($8 == "1|1") print $1,$2,$3,$4,$5,$6,"hap1,hap2",$8,$9}' results/integrate/mid.bed > results/integrate/mid2.bed
cat <(echo -e "#CHROM\tPOS\tID\tEND\tSVTYPE\tSVLEN\tHAP\tGT\tQUERY")  results/integrate/mid2.bed > results/integrate/LSG.bed

# each sv bed of hap1
$process_trans && grep 'TRANS' $hap1_sv_file > results/integrate/h1trans.bed
$process_sdr && grep 'SDR' $hap1_sv_file > results/integrate/h1sdr.bed
$process_dup && grep 'DUP' $hap1_sv_file > results/integrate/h1dup.bed
$process_highdup && grep 'HighDup' $hap1_sv_file > results/integrate/h1highdup.bed
$process_complex && grep 'COMPLEX' $hap1_sv_file > results/integrate/h1complex.bed

# each sv bed of hap2
if $process_trans; then
  grep 'TRANS' $hap2_sv_file | awk -F "\t" 'OFS="\t"{print $1,$2,$3,$4,$5,$6,$7,".|1",$9}' > results/integrate/h2trans.bed
fi
if $process_sdr; then
  grep 'SDR' $hap2_sv_file | awk -F "\t" 'OFS="\t"{print $1,$2,$3,$4,$5,$6,$7,".|1",$9}' > results/integrate/h2sdr.bed
fi
if $process_dup; then
  grep 'DUP' $hap2_sv_file | awk -F "\t" 'OFS="\t"{print $1,$2,$3,$4,$5,$6,$7,".|1",$9}' > results/integrate/h2dup.bed
fi
if $process_highdup; then
  grep 'HighDup' $hap2_sv_file | awk -F "\t" 'OFS="\t"{print $1,$2,$3,$4,$5,$6,$7,".|1",$9}' > results/integrate/h2highdup.bed
fi
if $process_complex; then
  grep 'COMPLEX' $hap2_sv_file | awk -F "\t" 'OFS="\t"{print $1,$2,$3,$4,$5,$6,$7,".|1",$9}' > results/integrate/h2complex.bed
fi

# merge sv bed 
sv_bed_files=""
$process_trans && sv_bed_files+="results/integrate/h1trans.bed results/integrate/h2trans.bed "
$process_sdr && sv_bed_files+="results/integrate/h1sdr.bed results/integrate/h2sdr.bed "
$process_dup && sv_bed_files+="results/integrate/h1dup.bed results/integrate/h2dup.bed "
$process_highdup && sv_bed_files+="results/integrate/h1highdup.bed results/integrate/h2highdup.bed "
$process_complex && sv_bed_files+="results/integrate/h1complex.bed results/integrate/h2complex.bed "

cat $sv_bed_files > results/integrate/sv.bed
\
awk '{print $0 "\t" NR}' results/integrate/sv.bed > results/integrate/sv.num.bed
bedtools coverage -a results/integrate/sv.num.bed -b results/integrate/h1_h2intersec.bed | awk '$14==1{print $0}' | awk 'OFS="\t"{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' | uniq > results/integrate/overlapped.bed
awk '{print $10}' results/integrate/overlapped.bed | uniq > results/integrate/overlap.num
end=$(wc -l < results/integrate/sv.bed)
seq 1 $end > results/integrate/all.num
sort results/integrate/all.num results/integrate/overlap.num | uniq -u | sort -n > results/integrate/sediff.num
paste <(awk 'OFS="\t"{print $1,$2,$3,$4,$5,$6,$7}' results/integrate/overlapped.bed) <(cut -f 8 results/integrate/overlapped.bed | sed 's/\./0/g') <(awk 'OFS="\t"{print $9}' results/integrate/overlapped.bed) > results/integrate/over.bed
awk 'NR==FNR {lines[$1]; next} FNR in lines' results/integrate/sediff.num results/integrate/sv.bed > results/integrate/sed.bed
cat results/integrate/sed.bed results/integrate/over.bed > results/integrate/LSGvarall.bed

cat results/integrate/LSG.bed results/integrate/LSGvarall.bed > $output_bed

rm -r results/integrate
