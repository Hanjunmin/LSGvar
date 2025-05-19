#!/bin/bash

# check params
if [ $# -lt 8 ]; then
  echo "Usage: $0 <pacigar.txt> <output_prefix> <input_vcf> <reference_fasta> <query_fasta> <SDR_output> <python_script> <bed_output> [variant_type]"
  echo "  variant_type: Optional - specify a single variant type to process (snv, ins, del, inv, trans, sdr, dup, highdup, complex)"
  echo "  If no variant_type is specified, all types will be processed."
  exit 1
fi

# get params
filein=$1       # cigar
fileout=$2      # file prefix
input_vcf=$3    # vcf file
reference_fasta=$4 # reference genome
query_fasta=$5    # query genome
SDR_output=$6     # sdr results
python_script=$7  # scripts
bed_output=$8     # bed file
variant_type=${9:-"all"}  # merge all variants

# VCF header
echo "Creating VCF header..."
rm -f vcf_header.txt
touch vcf_header.txt
echo -e "##fileformat=VCFv4.2" >>vcf_header.txt
time=$(date +'%Y%m%d')
echo -e "##fileDate=${time}" >>vcf_header.txt
echo -e "##source=LSGvar" >>vcf_header.txt
echo -e "##reference=file:$reference_fasta" >>vcf_header.txt
while IFS= read -r line; do
    contigname=$(echo "$line" |cut -f1)
    contiglen=$(echo "$line" |cut -f2)
    echo -e "##contig=<ID=${contigname},length=${contiglen}>" >>vcf_header.txt
done < "${reference_fasta}.fai"
echo -e "##INFO=<ID=ID,Number=1,Type=String,Description=\"Variant ID\">" >>vcf_header.txt
echo -e "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Variant type\">" >>vcf_header.txt
echo -e "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Variant length\">" >>vcf_header.txt
echo -e "##INFO=<ID=TIG_REGION,Number=.,Type=String,Description=\"Contig region where variant was found (one per alt with h1 before h2 for homozygous calls)\">" >>vcf_header.txt
echo -e "##INFO=<ID=QUERY_STRAND,Number=.,Type=String,Description=\"Strand of variant in the contig relative to the reference (order follows TIG_REGION)\">" >>vcf_header.txt
echo -e "##INFO=<ID=INNER_REF,Number=.,Type=String,Description=\"Inversion inner breakpoint in reference coordinates (order follows TIG_REGION)\">" >>vcf_header.txt
echo -e "##INFO=<ID=INNER_TIG,Number=.,Type=String,Description=\"Inversion inner breakpoint in contig coordinates (order follows TIG_REGION)\">" >>vcf_header.txt
echo -e "##INFO=<ID=HOM_REF,Number=.,Type=String,Description=\"Perfect breakpoint homology (SV sequence vs reference). Format 'X,Y' where X homology upstream, and Y is homology downstream. Homology vs reference is often better for DEL.\">" >>vcf_header.txt
echo -e "##INFO=<ID=HOM_TIG,Number=.,Type=String,Description=\"Perfect breakpoint homology (SV sequence vs contig). Format 'X,Y' where X homology upstream, and Y is homology downstream.  Homology vs contig is often better for INS.\">" >>vcf_header.txt
echo -e "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" >>vcf_header.txt
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample" >>vcf_header.txt

# Preprocess
echo "Preprocessing input file..."
paste $filein <(awk -F'[:-]' '{printf "%s\t%s\t%s\n", $2, $4, $1}' $filein) > addout.txt
awk -F'\t' 'NR>=2 {
    switch($5) {
        case "SNP_DEL":
              printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $11, $9 + $2 - 1, $9 + $2 + $4 - 1, $10 + $3, $10 + $3, $1, $4, $5, $6, $7, $8
            break
        case "SNP_INS":
                printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $11, $9 + $2 - 1, $9 + $2 - 1, $10 + $3, $10 + $3 + $4 - 1, $1, $4, $5, $6, $7, $8
            break
        case "INDEL_DEL":
            printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $11, $9 + $2 - 1, $9 + $2 + $4 - 1, $10 + $3, $10 + $3, $1, $4, $5, $6, $7, $8
            break
        case "INDEL_INS":
            printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $11, $9 + $2 - 1, $9 + $2 - 1, $10 + $3, $10 + $3 + $4 - 1, $1, $4, $5, $6, $7, $8
            break
        case "SNP":
            printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $11, $9 + $2, $9 + $2 + $4 - 1, $10 + $3, $10 + $3 + $4 - 1, $1, $4, $5, $6, $7, $8
            break
    }
}' addout.txt > CIGARend.txt

sed -i 's/SNP_DEL/DEL/g' CIGARend.txt
sed -i 's/INDEL_DEL/DEL/g' CIGARend.txt
sed -i 's/SNP_INS/INS/g' CIGARend.txt
sed -i 's/INDEL_INS/INS/g' CIGARend.txt

# Specific variants process
process_snv=false
process_ins=false
process_del=false
process_inv=false
process_trans=false
process_sdr=false
process_dup=false
process_highdup=false
process_complex=false

if [ "$variant_type" = "all" ]; then
  process_snv=true
  process_ins=true
  process_del=true
  process_inv=true
  process_trans=true
  process_sdr=true
  process_dup=true
  process_highdup=true
  process_complex=true
elif [ "$variant_type" = "snv" ]; then
  process_snv=true
elif [ "$variant_type" = "ins" ]; then
  process_ins=true
elif [ "$variant_type" = "del" ]; then
  process_del=true
elif [ "$variant_type" = "inv" ]; then
  process_inv=true
elif [ "$variant_type" = "trans" ]; then
  process_trans=true
elif [ "$variant_type" = "sdr" ]; then
  process_sdr=true
elif [ "$variant_type" = "dup" ]; then
  process_dup=true
elif [ "$variant_type" = "highdup" ]; then
  process_highdup=true
elif [ "$variant_type" = "complex" ]; then
  process_complex=true
else
  echo "Unknown variant type: $variant_type"
  echo "Valid types: snv, ins, del, inv, trans, sdr, dup, highdup, complex, all"
  exit 1
fi

# SNV
if $process_snv; then
  echo "Processing SNVs..."
  awk -F'\t' '$8=="SNP"{print$0}' CIGARend.txt >oursnv.txt
  file="oursnv.txt"
  paste ${file} <(awk '{print $1 "-" $2 "-" $8 "-" $9 "-" $10}' ${file}) <(awk -F'\t' '{print "ID=" $1 "-" $2 "-" $8 "-" $9 "-" $10 ";" "SVTYPE=" $8 ";" "TIG_REGION=" $1 ":" $4 "-" $5 ","  $1 ":" $4 "-" $5 ";" "QUERY_STRAND=" $11 ","$11}' ${file}) <(awk -F'\t' '{print "GT" }' ${file}) <(awk -F'\t' '{print "1|0" }' ${file}) >testbe.txt
  less testbe.txt|awk -F'\t' '{printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $12,$9,$10,".",".",$13,$14,$15}' >testchr1snvend.txt
  cat testchr1snvend.txt >endcigar.vcf
  sed 's/SNP/SNV/g' endcigar.vcf >end2cigar.vcf
  rm testchr1snvend.txt
fi

# DEL, INS, INV
if $process_del || $process_ins || $process_inv; then
  echo "Processing structural variants..."
  less -S $input_vcf | awk '$7 ~ "DEL" || $7 ~"INS" ||$7 ~"INV"{print $0}' | awk -F'\t' '($9!=0 || $10!=0){print $0}' >SDRend.txt
  sed -i 's\no\+\g' SDRend.txt
  sed -i 's/SDR_INV/INV/g' SDRend.txt
  sed -i 's/SV_INV/INV/g' SDRend.txt
  sed -i 's/SDR_INS/INS/g' SDRend.txt
  sed -i 's/SV_INS/INS/g' SDRend.txt
  sed -i 's/SDR_DEL/DEL/g' SDRend.txt
  sed -i 's/SV_DEL/DEL/g' SDRend.txt

  python $python_script --r $reference_fasta --q $query_fasta --i SDRend.txt --o $SDR_output
fi

# INS
if $process_ins; then
  echo "Processing insertions..."
  awk -F'\t' '$8=="INS"{print$0}' CIGARend.txt >ourins.txt
  file="ourins.txt"
  paste ${file} <(awk '{print $1 "-" $2 "-" $8 "-" $7}' ${file}) <(awk -F'\t' '{print "ID=" $1 "-" $2 "-" $8 "-" $7 ";" "SVTYPE=" $8 ";" "SVLEN=" $7 ";"  "TIG_REGION=" $1 ":" $4 "-" $5 ","  $1 ":" $4 "-" $5 ";" "QUERY_STRAND=" $11 ","$11}' ${file}) <(awk -F'\t' '{print "GT" }' ${file}) <(awk -F'\t' '{print "1|0" }' ${file}) >testbe.txt
  less testbe.txt|awk -F'\t' '{printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $12,$9,$10,".",".",$13,$14,$15}' >testchr1insend.txt
  rm ourins.txt
  rm testbe.txt
  cat <(less -S $SDR_output |awk 'index($3, "INS"){print$0}') testchr1insend.txt >ourinsend.txt
  rm testchr1insend.txt
fi

# DEL
if $process_del; then
  echo "Processing deletions..."
  awk -F'\t' '$8=="DEL"{print$0}' CIGARend.txt >ourdel.txt
  file="ourdel.txt"
  paste ${file} <(awk '{print $1 "-" $2 "-" $8 "-" $7}' ${file}) <(awk -F'\t' '{print "ID=" $1 "-" $2 "-" $8 "-" $7 ";" "SVTYPE=" $8 ";" "SVLEN=" "-" $7 ";"  "TIG_REGION=" $1 ":" $4 "-" $5 ","  $1 ":" $4 "-" $5 ";" "QUERY_STRAND=" $11 ","$11}' ${file}) <(awk -F'\t' '{print "GT" }' ${file}) <(awk -F'\t' '{print "1|0" }' ${file}) >testbe.txt
  less testbe.txt|awk -F'\t' '{printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $12,$9,$10,".",".",$13,$14,$15}' >testchr1delend.txt
  rm ourdel.txt
  rm testbe.txt
  cat <(less -S $SDR_output |awk 'index($3, "DEL"){print$0}') testchr1delend.txt >ourdelend.txt
  rm testchr1delend.txt
fi

# INV
if $process_inv; then
  echo "Processing inversions..."
  cat <(less -S $SDR_output |awk 'index($3, "INV"){print$0}') >ourinvend.txt
fi

# Merge
echo "Merging results..."
> variants.vcf
if $process_snv; then
  cat end2cigar.vcf >> variants.vcf
fi
if $process_ins; then
  cat ourinsend.txt >> variants.vcf
fi
if $process_del; then
  cat ourdelend.txt >> variants.vcf
fi
if $process_inv; then
  cat ourinvend.txt >> variants.vcf
fi

cat vcf_header.txt <(less -S variants.vcf) > $fileout
hap=$(basename "$filein" | cut -c 1-4)

# Bed results
echo "Creating BED files..."
> LSGvar.bed
echo -e "#CHROM\tPOS\tID\tEND\tSVTYPE\tSVLEN\tHAP\tGT\tQUERY" > LSGvar.bed

if $process_snv; then
  paste <(less -S variants.vcf |grep 'SNV' |awk -v hap=${hap} 'OFS="\t"{print $1,$2,$3,$2+1,"SNV",1,hap,"1|."}' |less -S) <(less -S variants.vcf |grep 'SNV' |awk -F'TIG_REGION=' '{print $2}' |less -S |awk -F',' '{print $1}') > SNV.bed
  cat SNV.bed >> LSGvar.bed
fi

if $process_ins; then
  paste <(less -S variants.vcf |grep 'INS' |awk -v hap=${hap} 'OFS="\t"{print $1,$2,$3}') <(less -S variants.vcf |grep 'INS' |awk -F'[-\t]' -v hap=${hap} 'OFS="\t"{print $2+1,"INS",$6,hap,"1|."}' |less -S) <(less -S variants.vcf |grep 'INS'|awk -F'TIG_REGION=' '{print $2}' |less -S |awk -F',' '{print $1}') > INS.bed
  cat INS.bed >> LSGvar.bed
fi

if $process_del; then
  paste <(less -S variants.vcf |grep 'DEL' |awk -v hap=${hap} 'OFS="\t"{print $1,$2,$3}') <(less -S variants.vcf |grep 'DEL' |awk -F'[-\t]' -v hap=${hap} 'OFS="\t"{print $2+$6,"DEL",$6,hap,"1|."}' |less -S) <(less -S variants.vcf|grep 'DEL' |awk -F'TIG_REGION=' '{print $2}' |less -S |awk -F',' '{print $1}') > DEL.bed
  cat DEL.bed >> LSGvar.bed
fi

if $process_trans; then
  less -S $input_vcf |grep 'TRANS' |awk -v hap=${hap} 'OFS="\t"{print $1,$2,$1"-"$2"-TRANS-"$9,$3,"TRANS",$9,hap,"1|.",$4":"$5"-"$6}' |less -S > TRANS.bed
  cat TRANS.bed >> LSGvar.bed
fi

if $process_sdr; then
  less -S $input_vcf |grep 'NM' |awk -v hap=${hap} 'OFS="\t"{print $1,$2,$1"-"$2"-SDR-"$9,$3,"SDR",$9,hap,"1|.",$4":"$5"-"$6}' |less -S > SDR.bed
  cat SDR.bed >> LSGvar.bed
fi

if $process_dup; then
  less -S $input_vcf |grep 'DUP' |awk -v hap=${hap} 'OFS="\t"{print $1,$2,$1"-"$2"-DUP-"$9,$3,"DUP",$9,hap,"1|.",$4":"$5"-"$6}' |less -S > DUP.bed
  cat DUP.bed >> LSGvar.bed
fi

if $process_highdup; then
  less -S $input_vcf |grep 'high-dup' |awk -v hap=${hap} 'OFS="\t"{print $1,$2,$1"-"$2"-HighDup",$3,"HighDup",".",hap,"1|.","."}' |less -S > Highdup.bed
  cat Highdup.bed >> LSGvar.bed
fi

if $process_inv; then
  less -S $input_vcf |grep 'INV' |awk -v hap=${hap} 'OFS="\t"{print $1,$2,$1"-"$2"-INV-"$9,$3,"INV",$9,hap,"1|.",$4":"$5"-"$6}' |less -S > INV.bed
  cat INV.bed >> LSGvar.bed
fi

if $process_complex; then
  less -S $input_vcf |grep 'COMPLEX' |awk -v hap=${hap} 'OFS="\t"{print $1,$2,$1"-"$2"-SDR_COMPLEX-"$9,$3,"SDR_COMPLEX",$9,hap,"1|.",$4":"$5"-"$6}' |less -S > complex.bed
  cat complex.bed >> LSGvar.bed
fi

less -S LSGvar.bed | awk 'OFS="\t"{print $1,$2,$4,$3,$5,$6,$7,$8,$9}' > $bed_output

rm -f SDRend.txt addout.txt *.bed *.vcf CIGARend.txt oursnv.txt ourinsend.txt ourdelend.txt ourinvend.txt

echo "Done! Processed variant type: $variant_type"
