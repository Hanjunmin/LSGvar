#!/bin/bash

if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <paf_file> <reverse_queries_file> <genome_path> <output_path>"
    exit 1
fi

INPUT="$1"
REVERSE_QUERIES="$2" 
GENOME_PATH="$3"     
OUTPUT_PATH="$4"     

awk '{
    if ($5 == "-") {
        sum[$1] += $10;
        len[$1] = $2;
    } else {
        len[$1] = $2;
    }
} 
END {
    for (q in len) {
        if (sum[q] > len[q]/2) {
            print q;
        }
    }
}' "$INPUT" > "$REVERSE_QUERIES"

RC_PRE_FASTA="rcpre.fasta"
RC_FASTA="rc.fasta"
OTHER_FASTA="other.fasta"

CHROMOSOMES=$(tr '\n' ',' < "$REVERSE_QUERIES" | sed 's/,$//')

seqkit grep -p "$CHROMOSOMES" "$GENOME_PATH" > "$RC_PRE_FASTA"
seqkit seq -r -p "$RC_PRE_FASTA" > "$RC_FASTA"
seqkit grep -v -p "$CHROMOSOMES" "$GENOME_PATH" > "$OTHER_FASTA"
cat "$RC_FASTA" "$OTHER_FASTA" > "$OUTPUT_PATH"
samtools faidx "$OUTPUT_PATH"

rm "$RC_PRE_FASTA" "$RC_FASTA" "$OTHER_FASTA"

echo "Finished! Reverse-complemented genome saved to: $OUTPUT_PATH"
