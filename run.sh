tool_path="/home/jmhan/SDR/"
ref_path="/home/jmhan/breakpoints/chm13v2.0.fa"
#hap1_path="/home/jmhan/SDR/genome/chimpanzee/mPanTro3.hap1.cur.20231122.fasta"
hap1_path="/home/jmhan/SDR/run_chimpanzee/genomedata/query_h1.fa"
#hap2_path="/home/jmhan/SDR/genome/chimpanzee/mPanTro3.hap2.cur.20231122.fasta"
hap2_path="/home/jmhan/SDR/run_chimpanzee/genomedata/query_h2.fa"
ref_species="human"
query_species=""
species=("human" "chimpanzee")
Dct=TRUE
nowdic=$(pwd)
mappingtsv="${nowdic}/chromosome_mapping.tsv"


nohup minimap2 -t 12 -cx asm20 --secondary=no --eqx -Y -K 8G -s 1000 $ref_path $hap1_path  -o hm_prihap1.paf &
nohup minimap2 -t 12 -cx asm20 --secondary=no --eqx -Y -K 8G -s 1000 $ref_path $hap2_path  -o hm_prihap2.paf &
python "${tool_path}scripts/one2multi_filter.py" -m ${mappingtsv} -f  hm_prihap1.paf -1 6 -2 1 |awk '{print $6,$7,$1,$2}' >p_c_chrlen1.txt
python "${tool_path}scripts/one2multi_filter.py" -m ${mappingtsv} -f  hm_prihap2.paf -1 6 -2 1 |awk '{print $6,$7,$1,$2}' >p_c_chrlen2.txt

mkdir saffireh1
chmod +x "${tool_path}scripts/one2multi_filter.py"
python "${tool_path}scripts/one2multi_filter.py" -m ${mappingtsv} -f  hm_prihap1.paf -1 6 -2 1 \
	| rustybam trim-paf \
	| rustybam break-paf --max-size 5000 \
	| rustybam filter --paired-len 100000 \
	| awk '$10 >= 20000' \
	> hm_prihap1.flt.paf

Rscript "${tool_path}scripts/chaos_filt.r"  "${nowdic}/p_c_chrlen1.txt" "${nowdic}/saffireh1/" "${nowdic}/hm_prihap1.flt.paf" "${nowdic}/afterchaos_hap1.flt.paf"  200000 200000
awk 'BEGIN{OFS="\t"}{print $6, $8, $9, ($8+$9)/2, $1, $3, $4, ($3+$4)/2, $5}' "${nowdic}/afterchaos_hap1.flt.paf"  |sort -k1,1 -k2,2n > "${nowdic}/h1syntenic_blocks.tsv"

cd ${nowdic}
mkdir saffireh2
chmod +x "${tool_path}scripts/one2multi_filter.py"
python "${tool_path}scripts/one2multi_filter.py" -m ${mappingtsv} -f  hm_prihap2.paf -1 6 -2 1 \
	| rustybam trim-paf \
	| rustybam break-paf --max-size 5000 \
	| rustybam filter --paired-len 100000 \
	| awk '$10 >= 20000' \
	> hm_prihap2.flt.paf

Rscript "${tool_path}scripts/chaos_filt.r"  "${nowdic}/p_c_chrlen2.txt" "${nowdic}/saffireh2/" "${nowdic}/hm_prihap2.flt.paf" "${nowdic}/afterchaos_hap2.flt.paf"  200000 200000
awk 'BEGIN{OFS="\t"}{print $6, $8, $9, ($8+$9)/2, $1, $3, $4, ($3+$4)/2, $5}' "${nowdic}/afterchaos_hap2.flt.paf"  |sort -k1,1 -k2,2n > "${nowdic}/h2syntenic_blocks.tsv"

python "${tool_path}scripts/CIGAR.py" --p 30 --r $ref_path --q $hap1_path --paf "${nowdic}/afterchaos_38_002ma.filt.paf" --o "${nowdic}/macigar.txt"
python "${tool_path}scripts/CIGAR.py" --p 30 --r $ref_path --q $hap2_path --paf "${nowdic}/afterchaos_38_002pa.filt.paf" --o "${nowdic}/pacigar.txt"


# SV
mkdir denSDRhap1
Rscript "${tool_path}scripts/denSDR.r"  "${tool_path}scripts/denSDRfun.r" "${nowdic}/masyntenic_blocks.tsv" "${nowdic}/p_c_chrlen1.txt" "${nowdic}/denSDRhap1/"
mkdir denSDRhap2
Rscript "${tool_path}scripts/denSDR.r"  "${tool_path}scripts/denSDRfun.r" "${nowdic}/pasyntenic_blocks.tsv" "${nowdic}/p_c_chrlen2.txt" "${nowdic}/denSDRhap2/"
cat ${nowdic}/denSDRhap2/*end.tsv > "${nowdic}/denSDRhap2/SDRall.txt"
cat ${nowdic}/denSDRhap1/*end.tsv > "${nowdic}/denSDRhap1/SDRall.txt"
## 2vcf
#bash "${tool_path}scripts/cigar2vcf.sh" "${nowdic}/macigar.txt" "${nowdic}/macigarsdr.vcf" "${nowdic}/denSDRhap1/SDRall.txt" $ref_path $hap1_path "${nowdic}/macigarsdr.txt" "${tool_path}scripts/SDR_vcf.py" "4.2.1"  ## 4.2.1 /0.6
#bash "${tool_path}scripts/cigar2vcf.sh" "${nowdic}/pacigar.txt" "${nowdic}/pacigarsdr.vcf" "${nowdic}/denSDRhap2/SDRall.txt" $ref_path $hap2_path "${nowdic}/pacigarsdr.txt" "${tool_path}scripts/SDR_vcf.py" "4.2.1"    ## 4.2.1 /0.6

# new38cigarma="/home/jmhan/SDR/HG002/integrate/our/vcfsplit/compare_addcen/dupcorrect/38/ma/38_002_macigiarout.txt"

# new38cigarpa="/home/jmhan/SDR/HG002/integrate/our/vcfsplit/compare_addcen/dupcorrect/38/pa/38_002_pacigiarout.txt"

# new37cigarma="/home/jmhan/SDR/HG002/integrate/our/vcfsplit/compare_addcen/dupcorrect/19/ma/37_002_macigiarout.txt"
# new37cigarpa="/home/jmhan/SDR/HG002/integrate/our/vcfsplit/compare_addcen/dupcorrect/19/pa/37_002_pacigiarout.txt"


bash "${tool_path}scripts/dup_filt.sh" "${nowdic}/afterchaos_38_002ma.filt.paf"  "${nowdic}/macigar.txt" "${nowdic}/macigarout.txt"  
bash "${tool_path}scripts/dup_filt.sh" "${nowdic}/afterchaos_38_002pa.filt.paf"  "${nowdic}/pacigar.txt" "${nowdic}/pacigarout.txt"  

bash "${tool_path}scripts/cigar2vcf.sh" "${new38cigarma}" "${nowdic}/macigarsdr.vcf" "${nowdic}/denSDRhap1/SDRall.txt" $ref_path $hap1_path "${nowdic}/macigarsdr.txt" "${tool_path}scripts/SDR_vcf.py" "0.6"  ## 4.2.1 /0.6
bash "${tool_path}scripts/cigar2vcf.sh" "${new38cigarpa}" "${nowdic}/pacigarsdr.vcf" "${nowdic}/denSDRhap2/SDRall.txt" $ref_path $hap2_path "${nowdic}/pacigarsdr.txt" "${tool_path}scripts/SDR_vcf.py" "0.6"    ## 4.2.1 /0.6


## phenotype
#split
bash "${tool_path}scripts/splitfile.sh" "${nowdic}/ma" "${nowdic}/macigarsdr.vcf" 
bash "${tool_path}scripts/splitfile.sh" "${nowdic}/pa" "${nowdic}/pacigarsdr.vcf"
#integrate
mkdir integrate && cd  integrate
bash "${tool_path}scripts/phenotype.sh"
