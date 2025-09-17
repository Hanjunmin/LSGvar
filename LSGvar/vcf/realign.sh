#!/bin/bash
usage() {
    echo "Usage: $0 -i input_file -r reference_fasta -q query_fasta -o output_file -t threads -c chunk_size -p paf_file"
    echo "Options:"
    echo "  -i  Input TSV file"
    echo "  -r  Reference genome FASTA"
    echo "  -q  Query haplotype FASTA"
    echo "  -o  Output file"
    echo "  -t  Threads"
    echo "  -c  Chunk_size"
    echo "  -p  paf_file"
    exit 1
}

while getopts "i:r:q:o:t:c:p:" opt; do
    case $opt in
        i) input="$OPTARG" ;;
        r) ref_path="$OPTARG" ;;
        q) hap_path="$OPTARG" ;;
        o) output="$OPTARG" ;;
        t) threads="$OPTARG" ;;
        c) chunk_size="$OPTARG" ;;
        p) paf_file="$OPTARG" ;;
        *) usage ;;
    esac
done

if [[ -z "$input" || -z "$ref_path" || -z "$hap_path" || -z "$output" || -z "$threads" || -z "$chunk_size" || -z "$paf_file" ]]; then
    echo "Error: Missing required parameters!"
    usage
fi

tmpdir=$(mktemp -d -p . -t tmp.XXXXXXXXXX)

max_threads=$threads
chunk_size=$chunk_size

process_chunk() {
    local chunk_file="$1"
    local chunk_tmpdir=$(mktemp -d -p "$tmpdir")
    local chunk_output="$chunk_tmpdir/result.tsv"
    
    while IFS=$'\t' read -r -a cols; do
        newline=("${cols[@]}")

        if [[ "${cols[0]}" == "ref_chr" ]]; then
            printf "%s\n" "$(IFS=$'\t'; echo "${newline[*]}")" >> "$chunk_output"
            continue
        fi

        # SDR_NM/SV_NM
        if [[ "${cols[6]}" == "SDR_NM" || "${cols[6]}" == "SV_NM" ]]; then
            ref_chr="${cols[0]}"
            ref_start="${cols[1]}"
            ref_end="${cols[2]}"
            q_chr="${cols[3]}"
            q_start="${cols[4]}"
            q_end="${cols[5]}"
            ref_len="${cols[8]}"
            q_len="${cols[9]}"

            len=$(( ref_len < q_len ? ref_len : q_len ))

            if (( len < 100 )); then
                printf "%s\n" "$(IFS=$'\t'; echo "${newline[*]}")" >> "$chunk_output"
                continue
            fi

            if (( len >= 5000000 )) || 
                (( ref_len > q_len && $(awk -v q="$q_len" -v r="$ref_len" 'BEGIN { print (q/r < 0.8) }') == 1 )) ||
                (( ref_len <= q_len && $(awk -v q="$q_len" -v r="$ref_len" 'BEGIN { print (r/q < 0.8) }') == 1 ))
            then
                printf "%s\n" "$(IFS=$'\t'; echo "${newline[*]}")" >> "$chunk_output"
                continue
            fi

            # inversion re-identification
            ref_fa="$chunk_tmpdir/ref_${ref_chr}_${ref_start}-${ref_end}.fa"
            query_fa="$chunk_tmpdir/query_${q_chr}_${q_start}-${q_end}.fa"
            
            if samtools faidx "$ref_path" "${ref_chr}:${ref_start}-${ref_end}" > "$ref_fa" 2>/dev/null &&
               samtools faidx "$hap_path" "${q_chr}:${q_start}-${q_end}" > "$query_fa" 2>/dev/null
            then
                minimap2 -t 24 -cx asm20 --eqx --secondary=no "$ref_fa" "$query_fa" > half/aln.paf 2>> "$chunk_tmpdir/minimap2.log"
                cat half/aln.paf >> half/minimap2.paf
                if [[ -s "half/aln.paf" ]]; then
                    total_matches=0
                    has_negative=false
                    has_positive=false
                    all_negative=true

                    while IFS=$'\t' read -ra paf_cols; do
                        if [[ "${paf_cols[4]}" == "-" ]]; then
                            has_negative=true
                        else
                            has_positive=true
                            all_negative=false
                        fi
                        ((total_matches += paf_cols[9]))  
                    done < "half/aln.paf"
    
                    if [[ "$has_negative" == "true" ]] && [[ "$all_negative" == "false" ]] && (( total_matches > len * 5 / 10 )); then
                        while IFS=$'\t' read -ra paf_cols; do
                            if [[ "${paf_cols[4]}" == "-" ]]; then
                                if [[ "${cols[6]}" == "SDR_NM" ]]; then
                                    newline[6]="SDR_NM(INV)"
                                elif [[ "${cols[6]}" == "SV_NM" ]]; then
                                    newline[6]="SV_NM(INV)"
                                fi
                                printf "%s\n" "$(IFS=$'\t'; echo "${newline[*]}")" >> "$chunk_output"
                                break  
                            fi
                        done < "half/aln.paf"
                    fi
    
                    if [[ "$all_negative" == "true" ]] && (( total_matches > len * 8 / 10 )); then
                        if [[ "${cols[6]}" == "SDR_NM" ]]; then
                            newline[6]="SDR_INV"
                        elif [[ "${cols[6]}" == "SV_NM" ]]; then
                            newline[6]="SV_INV"
                        fi
                        newline[7]="-"
                        printf "%s\n" "$(IFS=$'\t'; echo "${newline[*]}")" >> "$chunk_output"
                    fi
                fi
            else
                echo "ERROR: Failed to extract sequences for ${ref_chr}:${ref_start}-${ref_end}" >&2
            fi

        # SDR_INS 
        elif [[ "${cols[6]}" == "SDR_INS" ]]; then
            ref_chr="${cols[0]}"
            ref_start="${cols[1]}"
            local ins_result="$chunk_tmpdir/ins.result"

            awk -v chr="$ref_chr" '$6 == chr' "$paf_file" | cut -f6,8-9 | sort -k2,2n |
            bedtools merge -i - -d 500 |  
            bedtools closest -a <(echo -e "$ref_chr\t$ref_start\t$((ref_start+1))") -b - -D a |
            awk -v OFS="\t" '{
                dist_start = ($3 - $6) < 0 ? ($6 - $3) : ($3 - $6);
                dist_end = ($3 - $7) < 0 ? ($7 - $3) : ($3 - $7);
                min_dist = (dist_start < dist_end) ? dist_start : dist_end;
                
                if (min_dist <= 10) print "YES";
                else print "NO";
            }' > "$ins_result"
            
            if [[ -s "$ins_result" ]] && grep -q "YES" "$ins_result"; then
                  printf "%s\n" "$(IFS=$'\t'; echo "${newline[*]}")" >> "$chunk_output"
                  continue
            else
                  continue
            fi
        fi

        printf "%s\n" "$(IFS=$'\t'; echo "${newline[*]}")" >> "$chunk_output"
        
    done < "$chunk_file"
}

export -f process_chunk
export ref_path hap_path tmpdir paf_file
split -l $chunk_size --numeric-suffixes --additional-suffix=".tsv" "$input" "$tmpdir/chunk_"
find "$tmpdir" -name "chunk_*.tsv" | sort | xargs -P $max_threads -I {} bash -c 'process_chunk "{}"'
find "$tmpdir" -name "result.tsv" -exec cat {} + > "$output.tmp"
cat "$output.tmp" | uniq > "$output"
rm -r "$tmpdir" "$output.tmp"
