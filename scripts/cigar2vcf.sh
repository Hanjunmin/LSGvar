## header
rm vcf_header.txt
touch vcf_header.txt
echo -e "##fileformat=VCFv4.2" >>vcf_header.txt
time=$(date +'%Y%m%d')
echo -e "##fileDate=${time}" >>vcf_header.txt
echo -e "##source=LSGvar" >>vcf_header.txt
echo -e "##reference=file:$4" >>vcf_header.txt
while IFS= read -r line; do
    contigname=$(echo "$line" |cut -f1)
    contiglen=$(echo "$line" |cut -f2)
    echo -e "##contig=<ID=${contigname},length=${contiglen}>" >>vcf_header.txt
done < "${4}.fai"
echo -e "##INFO=<ID=ID,Number=1,Type=String,Description="Variant ID">" >>vcf_header.txt
echo -e "##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Variant type">" >>vcf_header.txt
echo -e "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Variant length">" >>vcf_header.txt
echo -e "##INFO=<ID=TIG_REGION,Number=.,Type=String,Description="Contig region where variant was found \(one per alt with h1 before h2 for homozygous calls\)">" >>vcf_header.txt
echo -e "##INFO=<ID=QUERY_STRAND,Number=.,Type=String,Description="Strand of variant in the contig relative to the reference \(order follows TIG_REGION\)">" >>vcf_header.txt
echo -e "##INFO=<ID=INNER_REF,Number=.,Type=String,Description="Inversion inner breakpoint in reference coordinates \(order follows TIG_REGION\)">" >>vcf_header.txt
echo -e "##INFO=<ID=INNER_TIG,Number=.,Type=String,Description="Inversion inner breakpoint in contig coordinates \(order follows TIG_REGION\)">" >>vcf_header.txt
echo -e "##INFO=<ID=HOM_REF,Number=.,Type=String,Description="Perfect breakpoint homology \(SV sequence vs reference\). Format 'X,Y' where X homology upstream, and Y is homology downstream. Homology vs reference is often better for DEL.">" >>vcf_header.txt
echo -e "##INFO=<ID=HOM_TIG,Number=.,Type=String,Description="Perfect breakpoint homology \(SV sequence vs contig\). Format 'X,Y' where X homology upstream, and Y is homology downstream.  Homology vs contig is often better for INS.">" >>vcf_header.txt
echo -e "##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">" >>vcf_header.txt
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample" >>vcf_header.txt
    
##
filein=$1 #pacigar.txt
fileout=$2
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
}' addout.txt  > CIGARend.txt

sed -i 's/SNP_DEL/DEL/g' CIGARend.txt
sed -i 's/INDEL_DEL/DEL/g' CIGARend.txt
sed -i 's/SNP_INS/INS/g' CIGARend.txt
sed -i 's/INDEL_INS/INS/g' CIGARend.txt


awk -F'\t' '$8=="SNP"{print$0}' CIGARend.txt >oursnv.txt
file="oursnv.txt"
paste ${file} <(awk '{print $1 "-" $2 "-" $8 "-" $9 "-" $10}' ${file}) <(awk -F'\t' '{print "ID=" $1 "-" $2 "-" $8 "-" $9 "-" $10 ";" "SVTYPE=" $8 ";" "TIG_REGION=" $1 ":" $4 "-" $5 ","  $1 ":" $4 "-" $5 ";" "QUERY_STRAND=" $11 ","$11}' ${file}) <(awk -F'\t' '{print "GT" }' ${file}) <(awk -F'\t' '{print "1|0" }' ${file}) >testbe.txt
#less -S chr1snv.vcf | cut -f 8 | sed 's/;/\t/g' |  cut -f 4 |  grep -n  "-" | less -S
less testbe.txt|awk -F'\t' '{printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $12,$9,$10,".",".",$13,$14,$15}' >testchr1snvend.txt
cat testchr1snvend.txt >endcigar.vcf
sed 's/SNP/SNV/g' endcigar.vcf >end2cigar.vcf



less -S $3 |awk '$7 ~ "DEL" || $7 ~"INS" ||$7 ~"INV"{print $0}' |awk -F'\t' '($9!=0 || $10!=0){print $0}' >SDRend.txt
sed -i 's\no\+\g' SDRend.txt
sed -i 's/SDR_INV/INV/g' SDRend.txt
sed -i 's/SV_INV/INV/g' SDRend.txt
sed -i 's/SDR_INS/INS/g' SDRend.txt
sed -i 's/SV_INS/INS/g' SDRend.txt
sed -i 's/SDR_DEL/DEL/g' SDRend.txt
sed -i 's/SV_DEL/DEL/g' SDRend.txt

python $7 --r $4 --q $5 --i SDRend.txt --o $6 ##生成了SDR.txt，这是我们的vcf，然后再整合
##INS
awk -F'\t' '$8=="INS"{print$0}' CIGARend.txt >ourins.txt
file="ourins.txt"
paste ${file} <(awk '{print $1 "-" $2 "-" $8 "-" $7}' ${file}) <(awk -F'\t' '{print "ID=" $1 "-" $2 "-" $8 "-" $7 ";" "SVTYPE=" $8 ";" "SVLEN=" $7 ";"  "TIG_REGION=" $1 ":" $4 "-" $5 ","  $1 ":" $4 "-" $5 ";" "QUERY_STRAND=" $11 ","$11}' ${file}) <(awk -F'\t' '{print "GT" }' ${file}) <(awk -F'\t' '{print "1|0" }' ${file}) >testbe.txt
less testbe.txt|awk -F'\t' '{printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $12,$9,$10,".",".",$13,$14,$15}' >testchr1insend.txt
rm ourins.txt
rm testbe.txt
cat <(less -S $6 |awk 'index($3, "INS"){print$0}') testchr1insend.txt >ourinsend.txt



awk -F'\t' '$8=="DEL"{print$0}' CIGARend.txt >ourdel.txt
file="ourdel.txt"
paste ${file} <(awk '{print $1 "-" $2 "-" $8 "-" $7}' ${file}) <(awk -F'\t' '{print "ID=" $1 "-" $2 "-" $8 "-" $7 ";" "SVTYPE=" $8 ";" "SVLEN=" "-" $7 ";"  "TIG_REGION=" $1 ":" $4 "-" $5 ","  $1 ":" $4 "-" $5 ";" "QUERY_STRAND=" $11 ","$11}' ${file}) <(awk -F'\t' '{print "GT" }' ${file}) <(awk -F'\t' '{print "1|0" }' ${file}) >testbe.txt
less testbe.txt|awk -F'\t' '{printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $12,$9,$10,".",".",$13,$14,$15}' >testchr1delend.txt
rm ourdel.txt
rm testbe.txt
cat <(less -S $6 |awk 'index($3, "DEL"){print$0}') testchr1delend.txt >ourdelend.txt





cat end2cigar.vcf  ourinsend.txt ourdelend.txt>hg002cigar.vcf
cat vcf_header.txt <(less -S hg002cigar.vcf)  >$fileout
hap=$(basename "${nowdic}/h1cigarout.txt"  |cut -c 1-2)

paste  <(less -S hg002cigar.vcf |grep 'SNV' |awk -v hap=${hap} 'OFS="\t"{print $1,$2,$3,$2+1,"SNV",1,hap,"1|."}' |less -S)    <(less -S hg002cigar.vcf |grep 'SNV' |awk -F'TIG_REGION=' '{print $2}' |less -S |awk -F',' '{print $1}') >SNV.bed
paste <(less -S hg002cigar.vcf |grep 'INS' |awk  -v hap=${hap} 'OFS="\t"{print $1,$2,$3}')  <(less -S hg002cigar.vcf |grep 'INS' |awk -F'[-\t]' -v hap=${hap} 'OFS="\t"{print $2+1,"INS",$6,hap,"1|."}' |less -S)    <(less -S hg002cigar.vcf |grep 'INS'|awk -F'TIG_REGION=' '{print $2}' |less -S |awk -F',' '{print $1}') >INS.bed
paste <(less -S hg002cigar.vcf |grep 'DEL' |awk  -v hap=${hap} 'OFS="\t"{print $1,$2,$3}') <(less -S hg002cigar.vcf |grep 'DEL' |awk -F'[-\t]' -v hap=${hap} 'OFS="\t"{print $2+$6,"DEL",$6,hap,"1|."}' |less -S)    <(less -S hg002cigar.vcf|grep 'DEL'  |awk -F'TIG_REGION=' '{print $2}' |less -S |awk -F',' '{print $1}') >DEL.bed
less -S $3  |grep 'TRANS' |awk -v hap=${hap} 'OFS="\t"{print $1,$2,$1"-"$2"-TRANS-"$9,$3,"TRANS",$9,hap,"1|.",$4":"$5"-"$6}' |less -S >TRANS.bed
less -S $3  |grep 'NM' |awk -v hap=${hap} 'OFS="\t"{print $1,$2,$1"-"$2"-SDR-"$9,$3,"SDR",$9,hap,"1|.",$4":"$5"-"$6}' |less -S >SDR.bed
less -S $3  |grep 'DUP' |awk -v hap=${hap} 'OFS="\t"{print $1,$2,$1"-"$2"-DUP-"$9,$3,"DUP",$9,hap,"1|.",$4":"$5"-"$6}' |less -S >DUP.bed
less -S $3  |grep 'high-dup' |awk -v hap=${hap} 'OFS="\t"{print $1,$2,$1"-"$2"-HighDup",$3,"HighDup",".",hap,"1|.","."}' |less -S >Highdup.bed
less -S $3  |grep 'INV' |awk -v hap=${hap} 'OFS="\t"{print $1,$2,$1"-"$2"-INV-"$9,$3,"INV",$9,hap,"1|.",$4":"$5"-"$6}' |less -S >INV.bed
less -S $3  |grep 'COMPLEX' |awk -v hap=${hap} 'OFS="\t"{print $1,$2,$1"-"$2"-COMPLEX-"$9,$3,"COMPLEX",$9,hap,"1|.",$4":"$5"-"$6}' |less -S >complex.bed
cat <(echo -e "#CHROM\tPOS\tID\tEND\tSVTYPE\tSVLEN\tHAP\tGT\tQUERY") SNV.bed INS.bed DEL.bed TRANS.bed  SDR.bed DUP.bed Highdup.bed INV.bed  complex.bed >LSGvar.bed
less -S LSGvar.bed |awk 'OFS="\t"{print $1,$2,$4,$3,$5,$6,$7,$8,$9}' > $8

rm SDRend.txt && rm addout.txt && rm *.bed && rm *.vcf && rm CIGARend.txt && rm oursnv.txt
# mark=$8
# if [ "$mark" = "4.2.1" ]; then
#     cat /home/jmhan/SDR/HG002/header.txt <(less -S hg002cigar.vcf |awk '$1!="chrX" && $1!="chrY"{print $0}' |awk 'NR>=3{print $0}')  >$fileout
# else
#     cat /home/jmhan/SDR/HG002/GRCH37/37header.txt <(less -S hg002cigar.vcf  |awk 'NR>=3{print $0}')  >$fileout
# fi


#/home/jmhan/SDR/HG002/header.txt    grch38的


## testchr1delend.txt  testchr1insend.txt end2cigar.vcf
