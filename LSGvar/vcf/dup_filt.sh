#!/bin/bash
paf=$1
file=$2
out=$3
script=$4

sort -k8,8n "$paf" |awk '{OFS=FS="\t"}{print $6,$8,$9}' > refpaf.region
awk '{count[$0]++} END {for (line in count) if (count[line] > 1) print line}' refpaf.region > duplicates.txt
bedtools intersect -a refpaf.region -b refpaf.region | awk 'NR==FNR{a[$0];next} !($0 in a)' refpaf.region - | sort -k1,1V -k2,2n -k3,3n | uniq | bedtools merge -i - > intersect.txt
cat duplicates.txt intersect.txt | sort -k1,1V -k2,2n -k3,3n > mergedup.txt

awk -F'[:-]' 'BEGIN {OFS="\t"}{print $0,$1,$2}' $file | awk '{OFS=FS="\t";print $0,$10+$2}' | nl > cigartestin.txt

awk -v OFS='\t' 'NR>1 {print $10, $11, $12, $1, $2}' cigartestin.txt > cigartestin.bed

bedtools intersect -a mergedup.txt -b cigartestin.bed -wa -wb > overlaps.txt
awk 'BEGIN {OFS="\t"} {if ($6 >= $2 && $6 <= $3) {print}}' overlaps.txt > overlaps_filt.txt

if [[ -s overlaps_filt.txt ]]; then
    python "$script" overlaps_filt.txt fp
    cat *_del.txt > del.txt
    
    if [ -s del.txt ]; then
        awk 'NR==FNR{d[$1]; next} !($1 in d)' del.txt cigartestin.txt | 
        cut -f 2-9 > "$out"
    else
        cut -f 2-9 cigartestin.txt > "$out"
    fi
else
    cut -f 2-9 cigartestin.txt > "$out"
fi

find . -maxdepth 1 -type f \( -name "refpaf.region" \
                             -o -name "mergedup.txt" \
                             -o -name "cigartestin.txt" \
                             -o -name "cigartestin.bed" \
                             -o -name "overlaps.txt" \
                             -o -name "overlaps_filt.txt" \
                             -o -name "del.txt" \
                             -o -name "duplicates.txt" \
                             -o -name "intersect.txt" \
                             -o -name "*_del.txt" \
                             -o -name "*_counts.txt" \) -delete
