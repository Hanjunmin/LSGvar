import subprocess
from Bio import SeqIO
import re
from multiprocessing import Pool
import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='tsv2vcf')
parser.add_argument('--r', help='Reference fasta')
parser.add_argument('--q', help='query fasta')
parser.add_argument('--i', help='SDR result file')
parser.add_argument('--o', help='output dir')
args = parser.parse_args()
QUERY =args.q
REF =args.r
paf=args.i



df = pd.read_csv(paf, sep='\t',header=None)
with open(args.o, 'w') as file:
	## INS
	ins_rows = df[df.iloc[:, 6] == 'INS']
	ins_rows.to_csv('ins.csv', sep='\t', header=None, index=False)
	subprocess.run(f"awk '{{print $4,$5-1,$6,\"\",\"\",$8}}' ins.csv|tr ' ' '\t' > query.bed", shell=True)
	subprocess.run(f"bedtools getfasta -fi {QUERY} -bed query.bed -s -fo que.fa", shell=True)
	subprocess.run(f"awk '{{print $1,$2-1,$3}}' ins.csv|tr ' ' '\t'  > ref.bed", shell=True)  ##如果参考基因组对应的是正链不进行反转，如果是负链就进行反转
	subprocess.run(f"bedtools getfasta -fi {REF} -bed ref.bed -fo ref.fa", shell=True)
	refsequences = {record.id: record.seq for record in SeqIO.parse('ref.fa', 'fasta')}
	querysequences = {record.id: record.seq for record in SeqIO.parse('que.fa', 'fasta')}
	for index, row in ins_rows.iterrows():
		ref_chr = row[0]
		ref_start = row[1]
		ref_end = row[2]
		que_chr = row[3]
		que_start = row[4]
		que_end = row[5]
		lenall = row[9]
		minus=row[7]
		selected_reffa = str(refsequences[ref_chr+":"+str(ref_start-1)+"-"+str(ref_end)])
		selected_quefa = str(querysequences[que_chr+":"+str(que_start-1)+"-"+str(que_end)+"("+minus+")"])
		three=ref_chr+"-"+str(ref_start)+"-"+"INS"+"-"+str(lenall)
		emp="."
		alllen=three+";"+"SVTYPE=INS;SVLEN="+str(lenall)+";TIG_REGION="+que_chr+":"+str(que_start)+"-"+str(que_end)+","+que_chr+":"+str(que_start)+"-"+str(que_end)+";"+"QUERY_STRAND=+,+"
		GT="GT"
		phase="1|0"
		file.write(f"{ref_chr}\t{ref_start}\t{three}\t{selected_reffa}\t{selected_quefa}\t{emp}\t{emp}\t{alllen}\t{GT}\t{phase}\n")  # 添加换行符，以便每行占一行
	## DEL
	del_rows = df[df.iloc[:, 6] == 'DEL']
	del_rows.to_csv('del.csv', sep='\t', header=None, index=False)
	subprocess.run(f"awk '{{print $4,$5-1,$6,\"\",\"\",$8}}' del.csv|tr ' ' '\t' > query.bed", shell=True)
	subprocess.run(f"bedtools getfasta -fi {QUERY} -bed query.bed -s -fo que.fa", shell=True)
	subprocess.run(f"awk '{{print $1,$2-1,$3}}' del.csv|tr ' ' '\t'  > ref.bed", shell=True)  ##如果参考基因组对应的是正链不进行反转，如果是负链就进行反转
	subprocess.run(f"bedtools getfasta -fi {REF} -bed ref.bed -fo ref.fa", shell=True)
	refsequences = {record.id: record.seq for record in SeqIO.parse('ref.fa', 'fasta')}
	querysequences = {record.id: record.seq for record in SeqIO.parse('que.fa', 'fasta')}
	for index, row in del_rows.iterrows():
		ref_chr = row[0]
		ref_start = row[1]
		ref_end = row[2]
		que_chr = row[3]
		que_start = row[4]
		que_end = row[5]
		lenall = row[8]
		minus=row[7]
		selected_reffa = str(refsequences[ref_chr+":"+str(ref_start-1)+"-"+str(ref_end)])
		selected_quefa = str(querysequences[que_chr+":"+str(que_start-1)+"-"+str(que_end)+"("+minus+")"])
		three=ref_chr+"-"+str(ref_start)+"-"+"DEL"+"-"+str(lenall)
		emp="."
		alllen=three+";"+"SVTYPE=DEL;SVLEN="+str(lenall)+";TIG_REGION="+que_chr+":"+str(que_start)+"-"+str(que_end)+","+que_chr+":"+str(que_start)+"-"+str(que_end)+";"+"QUERY_STRAND=+,+"
		GT="GT"
		phase="1|0"
		file.write(f"{ref_chr}\t{ref_start}\t{three}\t{selected_reffa}\t{selected_quefa}\t{emp}\t{emp}\t{alllen}\t{GT}\t{phase}\n")  # 添加换行符，以便每行占一行
