#set -e
paf=$1
file=$2
out=$3
# paf="/home/jmhan/SDR/HG002/afterchaos_38_002ma.filt.paf" #
# file="/home/jmhan/SDR/HG002/macigar.txt"  #
# out="38_002_macigiarout.txt" #


# paf="/home/jmhan/SDR/HG002/GRCH37/afterchaos_37_002pa.filt.paf" #
# file="//home/jmhan/SDR/HG002/GRCH37/pacigar.txt"  #
# out="37_002_pacigiarout.txt" #
less -S $paf |less -S |sort -k8,8n |less -S  |awk '{OFS=FS="\t"}{print $6,$8,$9}' |less -S >refpaf.region
awk 'NR==FNR{a[$0];next} !($0 in a)' refpaf.region <(bedtools intersect -a refpaf.region -b refpaf.region)  |sort -k1,1V -k2,2n -k3,3n |uniq|bedtools merge -i - |less -S >mergedup.txt
less -S $file | awk -F'[:-]' 'BEGIN {OFS="\t"}{print $0,$1,$2}' |less -S |awk '{OFS=FS="\t";print $0,$10+$2}' |nl >cigartestin.txt
cat /dev/null > del.txt


less -S mergedup.txt |less -S |while read -r line;do
echo $line
chr=$(echo $line |cut -d' ' -f 1)
start=$(echo $line |cut -d' ' -f 2)
end=$(echo $line |cut -d' ' -f 3)
less -S cigartestin.txt |awk -v chr="$chr" -v start="$start" -v end="$end" -F"\t" '$10==chr && $12>=start && $12<=end {print $0}' >test.txt
nr=$(less -S test.txt |wc -l)
if [ "$nr" -ne 0 ]; then
sampleuniq=$(less -S test.txt |cut -f 2 |sort |uniq |wc -l)
if [ "$sampleuniq" -ne 1 ];then
mapfile -t sampleuniqlist < <(less -S test.txt | cut -f 2 | sort | uniq)
for m in "${sampleuniqlist[@]}";do
seqstart=$(less -S test.txt |awk -v m="$m" '$2==m{print $0}'|sort -k12,12n |less -S |head -1 |cut -f 12)
seqend=$(less -S test.txt |awk -v m="$m" '$2==m{print $0}'|sort -k12,12n |less -S |tail -1 |cut -f 12)
seqnum=$(less -S test.txt |awk -v m="$m" '$2==m{print $0}' |wc -l)
if [ "$seqend" -ne "$seqstart" ];then
	freq=$(echo "scale=10; $seqnum / ($seqend - $seqstart)" | bc)
	echo $m,$freq
fi
done >sample_freq.txt
sam_len=$(less -S sample_freq.txt |wc -l)
if [ "$sam_len" -ne 0 ];then
mapfile -t sample_del < <(less -S sample_freq.txt | sort -t"," -k2,2g | awk -F"," 'NR>=2{print $1}')
for sample in "${sample_del[@]}";do
less -S test.txt |awk -v sample="$sample" '$2==sample{print $1}'
done >>del.txt
fi
fi
fi
done

if [ -s del.txt ]; then
    awk 'NR==FNR{seq[$1]; next} !($1 in seq)' del.txt cigartestin.txt |awk 'OFS=FS="\t"{print $2,$3,$4,$5,$6,$7,$8,$9}'  > $out
else
    cat cigartestin.txt |awk 'OFS=FS="\t"{print $2,$3,$4,$5,$6,$7,$8,$9}'  > $out
fi


rm cigartestin.txt && rm refpaf.region

## test
# less -S mergedup.txt | less -S | while read -r line; do
#     chr=$(echo $line | cut -d' ' -f 1)
#     start=$(echo $line | cut -d' ' -f 2)
#     end=$(echo $line | cut -d' ' -f 3)
#     less -S macigartestout.txt | awk -v chr="$chr" -v start="$start" -v end="$end" -F"\t" '$10==chr && $12>=start && $12<=end {print $0}' > test.txt
#     nr=$(less -S test.txt | wc -l)
#     if [ "$nr" -ne 0 ]; then
#         sampleuniq=$(less -S test.txt | cut -f 2 | sort | uniq | wc -l)
#         echo $line
#         echo $sampleuniq
#     fi
# done
