source config.json
nowdic=$(pwd)
## make new folder to store the genome
if [ ! -d "${nowdic}/genomedata" ]; 键，然后
	mkdir "${nowdic}/genomedata"
fi
##  1 minimap2
minimap2 -t 12 -cx asm20 --secondary=no --eqx -Y -K 8G -s 1000 $ref_path $hap1_path  -o init_1.paf
if [ ! -d "${nowdic}/saffireh1" ]; 键，然后
		mkdir "${nowdic}/saffireh1"
	fi
less init_1.paf|awk '{print $6,$7,$1,$2}' >p_c_chrleninit1.txt
Rscript "${tool_path}scripts/cen_filtonly.r"  "${nowdic}/p_c_chrleninit1.txt" "${nowdic}/saffireh1/" init_1.paf "${nowdic}/afterchaos_hap1init.flt.paf" ${centro} ${telome}
line_count=$(less "${nowdic}/add.csv" | wc -l)
if [ "$line_count" -ne 1 ]; 键，然后
    python "${tool_path}scripts/error_assem2fa.py" --f "$hap1_path" --i "${nowdic}/add.csv"
else
    echo "next"
fi
less "${nowdic}/minus.txt"|while read i; do samtools faidx ${hap1_path} $i >> ./genomedata/aa.fa; done
seqkit seq  -r -p ./genomedata/aa.fa >"./genomedata/h1.fa"
less "${nowdic}/plus.txt"|while read i; do samtools faidx ${hap1_path} $i >> "./genomedata/h1.fa"; done
if [ "$line_count" -ne 1 ]; 键，然后
    less "${nowdic}/add.csvoutput.fasta" >> "./genomedata/h1.fa"
else
    echo "next"
fi
rm ./genomedata/aa.fa
 
 ##  2 minimap2
if [ -n "$hap2_path" ]; 键，然后  ## two haplotypes
    minimap2 -t 12 -cx asm20 --secondary=no --eqx -Y -K 8G -s 1000 $ref_path $hap2_path  -o init_2.paf
    less init_2.paf|awk '{print $6,$7,$1,$2}' >p_c_chrleninit2.txt
	if [ ! -d "${nowdic}/saffireh2" ]; 键，然后
		mkdir "${nowdic}/saffireh2"
	fi
 less init_2.paf|awk '{print $6,$7,$1,$2}' >p_c_chrleninit2.txt
Rscript "${tool_path}scripts/cen_filtonly.r"  "${nowdic}/p_c_chrleninit2.txt" "${nowdic}/saffireh2/" init_2.paf "${nowdic}/afterchaos_hap2init.flt.paf" ${centro} ${telome}
line_count=$(less "${nowdic}/add.csv" | wc -l)
if [ "$line_count" -ne 1 ]; 键，然后
    python "${tool_path}scripts/error_assem2fa.py" --f "$hap2_path" --i "${nowdic}/add.csv"
else
    echo "next"
fi
less "${nowdic}/minus.txt"|while read i; do samtools faidx ${hap2_path} $i >> ./genomedata/aa.fa; done
seqkit seq  -r -p ./genomedata/aa.fa >"./genomedata/h2.fa"
less "${nowdic}/plus.txt"|while read i; do samtools faidx ${hap2_path} $i >> "./genomedata/h2.fa"; done
if [ "$line_count" -ne 1 ]; 键，然后
    less "${nowdic}/add.csvoutput.fasta" >> "./genomedata/h2.fa"
else
    echo "next"
fi
 
fi



## query反转



