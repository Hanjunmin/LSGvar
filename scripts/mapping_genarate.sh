rustybam trim-paf hm_prihap1.paf `#trims back alignments that align the same query sequence more than once` \
    | rustybam break-paf --max-size 5000 `#breaks the alignment into smaller pieces on indels of 5000 bases or more` \
    | rustybam orient `#orients each contig so that the majority of bases are forward aligned` \
    | rustybam filter --paired-len 100000 `#filters for query sequences that have at least 100,000 bases aligned to a target across all alignments.` \
    | rustybam stats --paf `#calculates statistics from the trimmed paf file` \
    > test.saffire

num=700
less -SH test.saffire | \
awk '$1 ~ /chr/ {print $1}' | \
sort | \
uniq | \
while read i; do
    awk -v val="$i" -v num="$num" '$1 == val {count[$6]++} END {output = val; for (key in count) if (count[key] > num) { sub(/[+-]/, "", key); output = output "\t" key } sub(/,$/, "", output); print output }' test.saffire | \
    awk -v OFS=',' '{if (NF > 1) {printf $1 "\t"; for (i=2; i<NF; i++) printf "%s,", $i; if (NF) printf "%s", $NF; print ""}}'; 
done \
>chromosome_mapping.tsv
