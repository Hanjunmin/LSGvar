file=$1
paste  <(less -S $file |grep 'SNV' |awk 'OFS="\t"{print $1,$2,$3,$2+1,"SNV",1,".",$10}' |less -S)    <(less -S $file |grep 'SNV' |awk -F'TIG_REGION=' '{print $2}' |less -S |awk -F',' '{print $1}') >SNV.bed
paste <(less -S $file |grep 'INS' |awk 'OFS="\t"{print $1,$2,$3}')  <(less -S $file |grep 'INS' |awk -F'[-\t]'  'OFS="\t"{print $2+1,"INS",$6,".",$18}' |less -S)    <(less -S $file |grep 'INS'|awk -F'TIG_REGION=' '{print $2}' |less -S |awk -F',' '{print $1}') >INS.bed
paste <(less -S $file |grep 'DEL' |awk  'OFS="\t"{print $1,$2,$3}') <(less -S $file |grep 'DEL' |awk -F'[-\t]'  'OFS="\t"{print $2+$6,"DEL",$6,".",$18}' |less -S)    <(less -S $file|grep 'DEL'  |awk -F'TIG_REGION=' '{print $2}' |less -S |awk -F',' '{print $1}') >DEL.bed
paste <(less -S $file |grep 'INV' |awk  'OFS="\t"{print $1,$2,$3}') <(less -S $file |grep 'INV' |awk -F'[-\t]' 'OFS="\t"{print $2+$6,"INV",$6,".",$18}' |less -S)    <(less -S $file|grep 'INV'  |awk -F'TIG_REGION=' '{print $2}' |less -S |awk -F',' '{print $1}') >INV.bed
cat SNV.bed INS.bed DEL.bed INV.bed >mid.bed
awk 'OFS="\t"{if ($8 == "0|1" || $8 == ".|1") print $1,$2,$3,$4,$5,$6,"h2",$8,$9; else if ($8 == "1|0" || $8 == "1|.") print $1,$2,$3,$4,$5,$6,"h1",$8,$9; else if ($8 == "1|1") print $1,$2,$3,$4,$5,$6,"h1,h2",$8,$9}' mid.bed >mid2.bed
cat <(echo -e "#CHROM\tPOS\tID\tEND\tSVTYPE\tSVLEN\tHAP\tGT\tQUERY")  mid2.bed >LSG.bed


cat $2 |grep 'TRANS' >h1trans.bed
cat $2 |grep 'SDR' >h1sdr.bed
cat $2 |grep 'DUP' >h1dup.bed
cat $2 |grep 'HighDup' >h1highdup.bed
cat $2 |grep 'COMPLEX' >h1complex.bed
cat $3 |grep 'TRANS' |awk -F "\t" 'OFS="\t"{print $1,$2,$3,$4,$5,$6,$7,".|1",$9}' >h2trans.bed
cat $3 |grep 'SDR' |awk -F "\t" 'OFS="\t"{print $1,$2,$3,$4,$5,$6,$7,".|1",$9}' >h2sdr.bed
cat $3 |grep 'DUP' |awk -F "\t" 'OFS="\t"{print $1,$2,$3,$4,$5,$6,$7,".|1",$9}'>h2dup.bed
cat $3 |grep 'HighDup'|awk -F "\t" 'OFS="\t"{print $1,$2,$3,$4,$5,$6,$7,".|1",$9}' >h2highdup.bed
cat $3 |grep 'COMPLEX' |awk -F "\t" 'OFS="\t"{print $1,$2,$3,$4,$5,$6,$7,".|1",$9}'>h2complex.bed
cat h1trans.bed h1sdr.bed h1dup.bed h1highdup.bed  h1complex.bed  h2trans.bed h2sdr.bed h2dup.bed h2highdup.bed h2complex.bed >sv.bed


awk '{print $0 "\t" NR}' sv.bed >sv.num.bed
bedtools coverage -a sv.num.bed -b h1_h2intersec.bed |less |awk '$14==1{print $0}' |awk 'OFS="\t"{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}'  |uniq > overlapped.bed
less -S overlapped.bed  |awk '{print $10}' |uniq >overlap.num
end=$(less sv.bed |wc -l)
seq 1 ${end} >all.num
sort all.num overlap.num| uniq -u |sort -n >sediff.num
paste <(less -S  overlapped.bed |awk 'OFS="\t"{print $1,$2,$3,$4,$5,$6,$7}') <(less -S  overlapped.bed |cut -f 8 |sed 's/\./0/g') <(less -S  overlapped.bed |awk 'OFS="\t"{print $9}')  >over.bed
awk 'NR==FNR {lines[$1]; next} FNR in lines' sediff.num sv.bed >sed.bed
cat sed.bed over.bed >LSGvarall.bed

cat  LSG.bed LSGvarall.bed >$4
