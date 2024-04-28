#hm_prihap1.flt为que ref
cd ${nowdic}
mkdir saffireh1
chmod +x "${tool_path}scripts/one2multi_filter.py"
python "${tool_path}scripts/one2multi_filter.py" -m ${mappingtsv} -f  hm_prihap1.paf -1 6 -2 1 \
	| rustybam trim-paf \
	| rustybam break-paf --max-size 5000 \
	| rustybam filter --paired-len 100000 \
	| rustybam orient \
	| awk '$10 >= 20000' \
	> hm_prihap1.flt.paf
if [ ! -d "${nowdic}/genomedata" ]; then
	mkdir "${nowdic}/genomedata"
fi
less -SH hm_prihap1.flt.paf | awk '{print $1}' | uniq | grep '-' | awk -F '-' '{print $1}' | sort | uniq > set1.txt
less -SH hm_prihap1.flt.paf | awk '{print $1}' | uniq | grep '+' | awk -F '+' '{print $1}' | sort | uniq > set2.txt

count=$(comm -12 set1.txt set2.txt | wc -l)
if [ $count -ne 0 ]; then
    awk 'BEGIN {OFS=FS="\t"} NR==FNR {a[$1"-"]; next} $1 in a {orig=$3; $1=substr($1, 1, length($1)-1)"+"; $3=$2-$4; $4=$2-orig; $5=($5=="+") ?"-" : "+";} 1' <(comm -12 set1.txt set2.txt) hm_prihap1.flt.paf > hm_prihap1.flta.paf
else
    less hm_prihap1.flt.paf > hm_prihap1.flta.paf
fi
less -SH hm_prihap1.flta.paf |awk '{print$1}' |uniq |grep '-'| awk -F '-' '{print$1}'|sort|uniq|less|while read i; do samtools faidx ${hap1_path} $i >> ./genomedata/aa.fa; done
seqkit seq  -r -p ./genomedata/aa.fa >./genomedata/query_h1.fa
less -SH hm_prihap1.flta.paf |awk '{print$1}' |uniq |grep '+'| awk -F '+' '{print$1}'|sort|uniq|less|while read i; do samtools faidx ${hap1_path} $i >> ./genomedata/query_h1.fa; done
## 由于原始序列可能大部分为负链，借用saffire中的orient参数得到大部分对应负链的query染色体，并将query原始序列反向互补
rm ./genomedata/aa.fa
#less hm_prihap1.flta.paf |awk 'gsub(/[-+]/,"",$1)'>Rinput_h1.paf


cd ${nowdic}
mkdir saffireh2
chmod +x "${tool_path}scripts/one2multi_filter.py"
python "${tool_path}scripts/one2multi_filter.py" -m ${mappingtsv} -f  hm_prihap2.paf -1 6 -2 1 \
	| rustybam trim-paf \
	| rustybam break-paf --max-size 5000 \
	| rustybam filter --paired-len 100000 \
	| rustybam orient \
	| awk '$10 >= 20000' \
	> hm_prihap2.flt.paf
less -SH hm_prihap2.flt.paf | awk '{print $1}' | uniq | grep '-' | awk -F '-' '{print $1}' | sort | uniq > set1.txt
less -SH hm_prihap2.flt.paf | awk '{print $1}' | uniq | grep '+' | awk -F '+' '{print $1}' | sort | uniq > set2.txt
count=$(comm -12 set1.txt set2.txt | wc -l)
if [ $count -ne 0 ]; then
    awk 'BEGIN {OFS=FS="\t"} NR==FNR {a[$1"-"]; next} $1 in a {orig=$3; $1=substr($1, 1, length($1)-1)"+"; $3=$2-$4; $4=$2-orig; $5=($5=="+") ?"-" : "+";} 1' <(comm -12 set1.txt set2.txt) hm_prihap2.flt.paf > hm_prihap2.flta.paf
else
    less hm_prihap2.flt.paf > hm_prihap2.flta.paf
fi
less -SH hm_prihap2.flta.paf |awk '{print$1}' |uniq |grep '-'| awk -F '-' '{print$1}'|sort|uniq|less|while read i; do samtools faidx ${hap2_path} $i >> ./genomedata/aa.fa; done
seqkit seq  -r -p ./genomedata/aa.fa >./genomedata/query_h2.fa
less -SH hm_prihap2.flta.paf |awk '{print$1}' |uniq |grep '+'| awk -F '+' '{print$1}'|sort|uniq|less|while read i; do samtools faidx ${hap2_path} $i >> ./genomedata/query_h2.fa; done
## 由于原始序列可能大部分为负链，借用saffire中的orient参数得到大部分对应负链的query染色体，并将query原始序列反向互补
rm ./genomedata/aa.fa
#less hm_prihap1.flta.paf |awk 'gsub(/[-+]/,"",$1)'>Rinput_h1.paf