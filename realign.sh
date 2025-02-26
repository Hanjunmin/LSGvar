#!/bin/bash

#Usage
usage() {
    echo "Usage: $0 -i input_file -r reference_fasta -q query_fasta -o output_file"
    exit 1
}

#params
while getopts "i:r:q:o:" opt; do
    case $opt in
        i) input="$OPTARG" ;;
        r) ref_path="$OPTARG" ;;
        q) hap_path="$OPTARG" ;;
        o) output="$OPTARG" ;;
        *) usage ;;
    esac
done

#recheck
if [[ -z "$input" || -z "$ref_path" || -z "$hap_path" || -z "$output" ]]; then
    echo "Error: no nessesory params!"
    usage
fi          

#temp directory
tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT

#deal with every row
while IFS=$'\t' read -r -a cols; do
    newline=("${cols[@]}")
    
    #extract the SDR_NM
    if [[ "${cols[6]}" == "SDR_NM" ]]; then
        #get ref and query coordinates
        ref_chr="${cols[0]}"
        ref_start="${cols[1]}"
        ref_end="${cols[2]}"
        q_chr="${cols[3]}"
        q_start="${cols[4]}"
        q_end="${cols[5]}"
        
        #extract the reference and query fasta
        ref_fa="$tmpdir/ref_${ref_chr}_${ref_start}_${ref_end}.fa"
        samtools faidx "$ref_path" "${ref_chr}:${ref_start}-${ref_end}" > "$ref_fa" 2>/dev/null
        
        query_fa="$tmpdir/query_${q_chr}_${q_start}_${q_end}.fa"
        samtools faidx "$hap_path" "${q_chr}:${q_start}-${q_end}" > "$query_fa" 2>/dev/null
        
        #redo the minimap alignment
        paf="$tmpdir/tmp.paf"
        minimap2 -t 12 -cx asm20 --secondary=no --eqx -Y -K 8G -s 1000 "$ref_fa" "$query_fa" -o "$paf" 2>/dev/null
        
        #check whether it is inversion
        if [[ -s "$paf" ]]; then
            echo "PAF file exists and is not empty."
            all_negative=true
            while IFS=$'\t' read -ra paf_cols; do
                strand="${paf_cols[4]}"
                if [[ "$strand" != "-" ]]; then
                    echo "Some alignments are on the positive strand."
                    all_negative=false
                    break
                fi
            done < "$paf"
    
            # If all alignments are negative strand, modify the row
            if [[ "$all_negative" == true ]]; then
                echo "All alignments are on the negative strand."
                newline[6]="SDR_INV"
                newline[7]="-"
            fi
        else
            echo "PAF file is empty or does not exist."
        fi
    fi
    
    (IFS=$'\t'; echo "${newline[*]}") >> "$output.tmp"
done < "$input"
mv "$output.tmp" "$output"

echo "Done!"
