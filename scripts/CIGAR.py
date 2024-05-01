## module load bedtools

import subprocess
from Bio import SeqIO
import re
from multiprocessing import Pool
import time
import argparse
parser = argparse.ArgumentParser(description='CIGAR-CALCULATE')
parser.add_argument('--p', help='number of process')
parser.add_argument('--r', help='Reference fasta')
parser.add_argument('--q', help='query fasta')
parser.add_argument('--paf', help='alignment paf file')
parser.add_argument('--o', help='output dir')
args = parser.parse_args()
print(args.p)


QUERY =args.q
REF =args.r
paf=args.paf
#paf="/home/jmhan/CIGAR/T2T/cigartest.paf"
subprocess.run(f"awk '{{print $1,$3,$4,\"\",\"\",$5}}' {paf}|tr ' ' '\t'  > query.bed", shell=True)
subprocess.run(f"bedtools getfasta -fi {QUERY} -bed query.bed -s -fo que.fa", shell=True)

subprocess.run(f"awk '{{print $6,$8,$9}}' {paf}|tr ' ' '\t'  > ref.bed", shell=True)  ##如果参考基因组对应的是正链不进行反转，如果是负链就进行反转
subprocess.run(f"bedtools getfasta -fi {REF} -bed ref.bed -fo ref.fa", shell=True)



refsequences = {record.id: record.seq for record in SeqIO.parse('ref.fa', 'fasta')}
querysequences = {record.id: record.seq for record in SeqIO.parse('que.fa', 'fasta')}


def process_line(line):
	columns = line.split('\t')
	minus = columns[4]  ## 判断正负链
	ref_chr = columns[5]
	ref_start = columns[7]
	print(ref_start)
	ref_end = columns[8]
	query_chr = columns[0]
	query_start = columns[2]
	query_end = columns[3]
	for xx in columns:
		if xx.startswith('cg'):
			CIGARstr = xx
			break
	index_of_first_digit = next((index for index, char in enumerate(CIGARstr) if char.isdigit()), None)
	CIGARstr =CIGARstr[index_of_first_digit:]
	selected_reffa = refsequences[ref_chr+":"+ref_start+"-"+ref_end]
	selected_quefa = querysequences[query_chr+":"+query_start+"-"+query_end+"("+minus+")"]
	id=ref_chr+":"+ref_start+"-"+ref_end+'_'+query_chr+":"+query_start+"-"+query_end
	matches = re.finditer(r'(\d+[DXI])', CIGARstr) # SNP
	alllen=int(query_end)-int(query_start)
	result_lines = []
	for match in matches:
		ref_posf=sum(list(map(int, re.findall(r'\d+',''.join(re.findall(r'(\d+[=XD])', CIGARstr[:match.start()]))))))
		query_posf=sum(list(map(int, re.findall(r'\d+',''.join(re.findall(r'(\d+[=XI])', CIGARstr[:match.start()]))))))
		if(match.group(1)=="1D"):
			ref_strseq=selected_reffa[ref_posf-1:ref_posf+1] ## ref_posf-1的位置是一样的
			que_strseq=selected_quefa[query_posf-1:query_posf] ##原来为空
			if(minus=="-"):
				result_lines.append(f"{id}\t{ref_posf+1}\t{alllen-query_posf}\t1\tSNP_DEL\t{ref_strseq}\t{que_strseq}\t{minus}\n") ## ref_posf ref_posf+1 query_posf query_posf 
			else:
				result_lines.append(f"{id}\t{ref_posf+1}\t{query_posf+1}\t1\tSNP_DEL\t{ref_strseq}\t{que_strseq}\t{minus}\n") ## ref_posf ref_posf+1 query_posf query_posf 
			
		if(match.group(1)=="1I"):
			ref_strseq=selected_reffa[ref_posf-1:ref_posf] #原先为空
			que_strseq=selected_quefa[query_posf-1:query_posf+1]
			
			if(minus=="-"):
				result_lines.append(f"{id}\t{ref_posf+1}\t{alllen-query_posf}\t1\tSNP_INS\t{ref_strseq}\t{que_strseq}\t{minus}\n")## ref_posf ref_posf query_posf query_posf+1 
			else:
				result_lines.append(f"{id}\t{ref_posf+1}\t{query_posf+1}\t1\tSNP_INS\t{ref_strseq}\t{que_strseq}\t{minus}\n")## ref_posf ref_posf query_posf query_posf+1 
			
		if(match.group(1)=="1X"):
			ref_strseq=selected_reffa[ref_posf:ref_posf+1]  ### 左闭右开，所以ref_posf个相等（最后一个位置是ref_posf-1），那新的变异应该是ref_posf,但是实际的位置应该是ref_posf+1
			que_strseq=selected_quefa[query_posf:query_posf+1]
			if(minus=="-"):
				result_lines.append(f"{id}\t{ref_posf+1}\t{alllen-query_posf}\t1\tSNP\t{ref_strseq}\t{que_strseq}\t{minus}\n") ## ref_posf ref_posf+1 query_posf query_posf+1 
			else:
				result_lines.append(f"{id}\t{ref_posf+1}\t{query_posf+1}\t1\tSNP\t{ref_strseq}\t{que_strseq}\t{minus}\n") ## ref_posf ref_posf+1 query_posf query_posf+1 
			
			
	# indel
	# ins
	matches = re.finditer(r'((?!1[DI])\d+[I])', CIGARstr)
	for match in matches:
		len=int(re.findall(r'\d+',match.group(1))[0])
		ref_posf=sum(list(map(int, re.findall(r'\d+',''.join(re.findall(r'(\d+[=XD])', CIGARstr[:match.start()]))))))
		query_posf=sum(list(map(int, re.findall(r'\d+',''.join(re.findall(r'(\d+[=XI])', CIGARstr[:match.start()]))))))
		ref_strseq=selected_reffa[ref_posf-1:ref_posf] ##原来为空
		que_strseq=selected_quefa[query_posf-1:query_posf+len]
		if(minus=="-"):
			result_lines.append(f"{id}\t{ref_posf+1}\t{alllen-query_posf-len+1}\t{len}\tINDEL_INS\t{ref_strseq}\t{que_strseq}\t{minus}\n")
		else:
			result_lines.append(f"{id}\t{ref_posf+1}\t{query_posf+1}\t{len}\tINDEL_INS\t{ref_strseq}\t{que_strseq}\t{minus}\n")
			
		
	matches = re.finditer(r'((?!1[DI])\d+[D])', CIGARstr)  ## del
	for match in matches:
		len=int(re.findall(r'\d+',match.group(1))[0])
		ref_posf=sum(list(map(int, re.findall(r'\d+',''.join(re.findall(r'(\d+[=XD])', CIGARstr[:match.start()]))))))
		query_posf=sum(list(map(int, re.findall(r'\d+',''.join(re.findall(r'(\d+[=XI])', CIGARstr[:match.start()]))))))
		ref_strseq=selected_reffa[ref_posf-1:ref_posf+len]
		que_strseq=selected_quefa[query_posf-1:query_posf] #原来为空
		if(minus=="-"):
			result_lines.append(f"{id}\t{ref_posf+1}\t{alllen-query_posf}\t{len}\tINDEL_DEL\t{ref_strseq}\t{que_strseq}\t{minus}\n")
		else:
			result_lines.append(f"{id}\t{ref_posf+1}\t{query_posf+1}\t{len}\tINDEL_DEL\t{ref_strseq}\t{que_strseq}\t{minus}\n")		
	

	matches = re.finditer(r'((?!1[DIX])\d+[X])', CIGARstr)
	for match in matches:
		len=int(re.findall(r'\d+',match.group(1))[0])
		j=0
		for i in range(len):
			ref_posf=sum(list(map(int, re.findall(r'\d+',''.join(re.findall(r'(\d+[=XD])', CIGARstr[:match.start()]))))))+j
			query_posf=sum(list(map(int, re.findall(r'\d+',''.join(re.findall(r'(\d+[=XI])', CIGARstr[:match.start()]))))))+j
			ref_strseq=selected_reffa[ref_posf:ref_posf+1] 
			que_strseq=selected_quefa[query_posf:query_posf+1]
			if(minus=="-"):
				result_lines.append(f"{id}\t{ref_posf+1}\t{alllen-query_posf}\t1\tSNP\t{ref_strseq}\t{que_strseq}\t{minus}\n")
			else:
				result_lines.append(f"{id}\t{ref_posf+1}\t{query_posf+1}\t1\tSNP\t{ref_strseq}\t{que_strseq}\t{minus}\n")		
			j=j+1

	with open("temp/"+ref_chr+":"+ref_start+"-"+ref_end+"_"+query_chr+":"+query_start+"-"+query_end+".cigar", 'w') as outfile:
		for line in result_lines:
			outfile.write(line)
	print("yes")
	return result_lines	

	

a=time.time()
with open(paf, 'r') as infile, open(args.o, 'w') as outfile:
	write=outfile.write(f"id\tref_start\tquery_start\tlen\ttype\tref_seq\tquery_seq\tQUERYSTRAND\n")
	with Pool(processes=int(args.p)) as pool:
		## 之前是30
		results = pool.map(process_line, infile)
	b=time.time()
	print((a-b)/60)


subprocess.run(f"rm query.bed", shell=True)
subprocess.run(f"rm que.fa", shell=True)

subprocess.run(f"rm ref.bed", shell=True)  ##如果参考基因组对应的是正链不进行反转，如果是负链就进行反转
subprocess.run(f"rm ref.fa", shell=True)
# with open(paf, 'r') as infile:
# 	all_lines = infile.readlines()



# for match in matches:
# 	print(match)
