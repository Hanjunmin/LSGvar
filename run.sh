source config.json
nowdic=$(pwd)

##-------------------------------------------------------- step1.minimap---------------------------------------------------------
minimap2 -t 12 -cx asm20 --secondary=no --eqx -Y -K 8G -s 1000 $ref_path $hap1_path  -o hm_prihap1.paf
python "${tool_path}/scripts/one2multi_filter.py" -m ${mappingtsvh1} -f  hm_prihap1.paf -1 6 -2 1 |awk '{print $6,$7,$1,$2}' >p_c_chrlen1.txt
python "${tool_path}/scripts/one2multi_filter.py" -m ${mappingtsvh1} -f  hm_prihap1.paf -1 6 -2 1 > hm_prihap1.flt.paf


##---------------------------------------------------------- step2.Filter----------------------------------------------------------
if [ ! -d "${nowdic}/saffireh1" ]; then
	mkdir "${nowdic}/saffireh1"
fi
if [ "$3" = "ctn" ]; then
    Rscript "${tool_path}/scripts/chaos_filt.r" "${nowdic}/p_c_chrlen1.txt" "${nowdic}/saffireh1/" "${nowdic}/hm_prihap1.flt.paf" "${nowdic}/afterchaos_hap1.flt.paf" $1 $2 $3
elif [ "$3" = "cts" ]; then
    Rscript "${tool_path}/scripts/chaos_filt.r" "${nowdic}/p_c_chrlen1.txt" "${nowdic}/saffireh1/" "${nowdic}/hm_prihap1.flt.paf" "${nowdic}/afterchaos_hap1.flt.paf" $1 $2 $3 ${centro} ${telome}
else
    echo "please input ctn/cts parameter"
fi
awk 'BEGIN{OFS="\t"}{print $6, $8, $9, ($8+$9)/2, $1, $3, $4, ($3+$4)/2, $5}' "${nowdic}/afterchaos_hap1.flt.paf"  |sort -k1,1 -k2,2n > "${nowdic}/h1syntenic_blocks.tsv"

##------------------------------------------------------------ step3.Another haplotype?----------------------------------------------------
if [ -n "$hap2_path" ]; then  ## two haplotypes
    minimap2 -t 12 -cx asm20 --secondary=no --eqx -Y -K 8G -s 1000 $ref_path $hap2_path  -o hm_prihap2.paf
	python "${tool_path}/scripts/one2multi_filter.py" -m ${mappingtsvh2} -f  hm_prihap2.paf -1 6 -2 1 |awk '{print $6,$7,$1,$2}' >p_c_chrlen2.txt
	python "${tool_path}/scripts/one2multi_filter.py" -m ${mappingtsvh2} -f  hm_prihap2.paf -1 6 -2 1  >hm_prihap2.flt.paf
	if [ ! -d "${nowdic}/saffireh2" ]; then
		mkdir "${nowdic}/saffireh2"
	fi
	if [ "$3" = "ctn" ]; then
	    Rscript "${tool_path}/scripts/chaos_filt.r" "${nowdic}/p_c_chrlen2.txt" "${nowdic}/saffireh2/" "${nowdic}/hm_prihap2.flt.paf" "${nowdic}/afterchaos_hap2.flt.paf" $1 $2 $3
	elif [ "$3" = "cts" ]; then
	    Rscript "${tool_path}/scripts/chaos_filt.r" "${nowdic}/p_c_chrlen2.txt" "${nowdic}/saffireh2/" "${nowdic}/hm_prihap2.flt.paf" "${nowdic}/afterchaos_hap2.flt.paf" $1 $2 $3 ${centro} ${telome}
	else
	    echo "please input ctn/cts parameter"
	fi

awk 'BEGIN{OFS="\t"}{print $6, $8, $9, ($8+$9)/2, $1, $3, $4, ($3+$4)/2, $5}' "${nowdic}/afterchaos_hap2.flt.paf"  |sort -k1,1 -k2,2n > "${nowdic}/h2syntenic_blocks.tsv"
fi


##------------------------------------------------------------------ step3.Cluster and call SV------------------------------------------------
[ -d denSDRhap1 ] || mkdir denSDRhap1
Rscript "${tool_path}/scripts/denSDR.r"  "${tool_path}/scripts/denSDRfun.r" "${nowdic}/h1syntenic_blocks.tsv" "${nowdic}/p_c_chrlen1.txt" "${nowdic}/denSDRhap1/"
if [ $? -eq 0 ]; then
    echo "hap1 CLUSTER Success" && rm "${nowdic}/h1syntenic_blocks.tsv" && rm "${nowdic}/p_c_chrlen1.txt"
else
    echo "hap1 CLUSTER Failure"
fi
cat ${nowdic}/denSDRhap1/*end.tsv >} "${nowdic/denSDRhap1/SDRall.txt"
bash "${tool_path}/realign.sh" -i "${nowdic/denSDRhap1/SDRall.txt" -r $ref_path -q $hap1_path -o "${nowdic/denSDRhap1/SDRall_final.txt"

if [ -n "$hap2_path" ]; then  ## two haplotypes
[ -d denSDRhap2 ] || mkdir denSDRhap2
Rscript "${tool_path}/scripts/denSDR.r"  "${tool_path}/scripts/denSDRfun.r" "${nowdic}/h2syntenic_blocks.tsv" "${nowdic}/p_c_chrlen2.txt" "${nowdic}/denSDRhap2/"
if [ $? -eq 0 ]; then
    echo "hap2 CLUSTER Success" && rm "${nowdic}/h2syntenic_blocks.tsv" && rm "${nowdic}/p_c_chrlen2.txt"
else
    echo "hap2 CLUSTER Failure"
fi
cat ${nowdic}/denSDRhap2/*end.tsv > "${nowdic}/denSDRhap2/SDRall.txt"
bash "${tool_path}/realign.sh" -i "${nowdic/denSDRhap2/SDRall.txt" -r $ref_path -q $hap2_path -o "${nowdic/denSDRhap2/SDRall_final.txt"
fi

##-------------------------------------------------------------------------- step4.CIGAR--------------------------------------------------------
[ -d temp ] || mkdir temp
python "${tool_path}/scripts/CIGAR.py" --r $ref_path --q $hap1_path --paf "${nowdic}/afterchaos_hap1.flt.paf" --o "${nowdic}/h1cigar.txt"
cat "${nowdic}/h1cigar.txt" temp/*.cigar >"${nowdic}/h1cigarend.txt"
rm temp/*.cigar
if [ -n "$hap2_path" ]; then  ## two haplotypes
	python "${tool_path}/scripts/CIGAR.py" --r $ref_path --q $hap2_path --paf "${nowdic}/afterchaos_hap2.flt.paf" --o "${nowdic}/h2cigar.txt"
	cat  "${nowdic}/h2cigar.txt" temp/*.cigar >"${nowdic}/h2cigarend.txt"
	rm temp/*.cigar
fi
[ -f "${nowdic}/ref.fa" ] && rm "${nowdic}/ref.fa"
[ -f "${nowdic}/que.fa" ] && rm "${nowdic}/que.fa"
[ -f "${nowdic}/query.bed" ] && rm "${nowdic}/query.bed"
[ -f "${nowdic}/ref.bed" ] && rm "${nowdic}/ref.bed"
##--------------------------------------------------------------------------- step5.dup_filt--------------------------------------------------
bash "${tool_path}/scripts/dup_filt.sh" "${nowdic}/afterchaos_hap1.flt.paf"  "${nowdic}/h1cigarend.txt" "${nowdic}/h1cigarout.txt"  
if [ -n "$hap2_path" ]; then  ## two haplotypes
	bash "${tool_path}/scripts/dup_filt.sh" "${nowdic}/afterchaos_hap2.flt.paf"  "${nowdic}/h2cigarend.txt" "${nowdic}/h2cigarout.txt"  
fi

## ------------------------------------------------------------------------step6.SV_INDEL_SNV2vcf---------------------------------------------
if [ ! -d "${nowdic}/results" ]; then
	mkdir "${nowdic}/results"
fi
bash "${tool_path}/scripts/cigar2vcf.sh" "${nowdic}/h1cigarout.txt" "${nowdic}/results/h1cigarsdr.vcf" "${nowdic}/denSDRhap1/SDRall_final.txt" $ref_path $hap1_path "${nowdic}/h1cigarsdr.txt" "${tool_path}/scripts/SDR_vcf.py" "${nowdic}/results/LSGvarend1.bed"
if [ -n "$hap2_path" ]; then  ## two haplotypes
bash "${tool_path}/scripts/cigar2vcf.sh" "${nowdic}/h2cigarout.txt" "${nowdic}/results/h2cigarsdr.vcf" "${nowdic}/denSDRhap2/SDRall_final.txt" $ref_path $hap2_path "${nowdic}/h2cigarsdr.txt" "${tool_path}/scripts/SDR_vcf.py" "${nowdic}/results/LSGvarend2.bed"
fi
[ -f "${nowdic}/inv.csv" ] && rm "${nowdic}/inv.csv"
[ -f "${nowdic}/ins.csv" ] && rm "${nowdic}/ins.csv"
[ -f "${nowdic}/del.csv" ] && rm "${nowdic}/del.csv"
##------------------------------------------------------------------------- step7.split and integrate------------------------------------------
bash "${tool_path}/scripts/splitfile.sh" "${nowdic}/h1" "${nowdic}/results/h1cigarsdr.vcf" 
if [ -n "$hap2_path" ]; then  ## two haplotypes
	bash "${tool_path}/scripts/splitfile.sh" "${nowdic}/h2" "${nowdic}/results/h2cigarsdr.vcf"
fi

# ## ----------------------------------------------------------------------------------------------------------
if [ -n "$hap2_path" ]; then
	[ -d results/integrate ] || mkdir -p results/integrate 
	cd  results/integrate
 	less -S "${nowdic}/afterchaos_hap1.flt.paf" |awk 'OFS="\t"{print $6,$8,$9}'|bedtools sort |bedtools merge -i - >h1_paf.bed
	less -S "${nowdic}/afterchaos_hap2.flt.paf" |awk 'OFS="\t"{print $6,$8,$9}'|bedtools sort |bedtools merge -i - >h2_paf.bed
	bedtools intersect -a h1_paf.bed -b h2_paf.bed  |bedtools sort >h1_h2intersec.bed
	bash "${tool_path}/scripts/phenotype.sh" "${nowdic}/h1" "${nowdic}/h2" ${ref_path} 
	bash "${tool_path}/scripts/vcf2bedGT.sh" "${nowdic}/results/sortLSGvarall.vcf.gz" "${nowdic}/results/LSGvarend1.bed" "${nowdic}/results/LSGvarend2.bed" "${nowdic}/results/LSGvar.bed"
fi

[ -f "${nowdic}/h2cigar.txt" ] && rm "${nowdic}/h2cigar.txt"
[ -f "${nowdic}/h1cigar.txt" ] && rm "${nowdic}/h1cigar.txt"
