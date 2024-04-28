import subprocess
from Bio import SeqIO
import re
from multiprocessing import Pool
import time
import argparse
parser = argparse.ArgumentParser(description='fasta')
parser.add_argument('--f', help='contig input')
parser.add_argument('--i', help='contig input')
args = parser.parse_args()
print(args.i)
refsequences = {record.id: record.seq for record in SeqIO.parse(args.f, 'fasta')}
import pandas as pd
df = pd.read_csv(args.i) 
outseq={}
for value in df["contig"].unique():
	selected_rows=df[df["contig"]==value]
	for index,row in selected_rows.iterrows():
		if row['strand']=="+":
			outseq[value+str(index)]=refsequences[value][row['start']:row['end']]
		else:
			outseq[value+str(index)]=refsequences[value][row['start']:row['end']].reverse_complement()

with open(args.i+"output.fasta", "w") as output_file:
	for name, seq in outseq.items():
		output_file.write(">" + name + "\n")
		output_file.write(str(seq) + "\n")