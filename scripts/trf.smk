import os
import sys
import pandas as pd
TRF_PATH = config.get("trf_path", "trf")
SAMTOOLS_PATH = config.get("samtools_path", "samtools")
ASSEMBLY = config["assembly"]
FAI = ASSEMBLY + ".fai"
PREFIX = config.get("prefix", "asm")
if not os.path.isfile(FAI):
    os.system(f"{SAMTOOLS_PATH} faidx {ASSEMBLY}")
if config["intervals"] and os.path.isfile(config["intervals"]):
    INTERVALS = [line.strip() for line in open(config["intervals"])]
else:
    INTERVALS = [line.strip() for line in os.popen(f"cut -f 1 {FAI}")]
wildcard_constraints:
	interval = "|".join(INTERVALS)

rule all:
	input:
		f"{PREFIX}telo.txt"
rule split_FASTA:
	input:
		ASSEMBLY
	output:
		temp("temp/{interval}.fasta")
	threads:
		1
	shell:"""
{SAMTOOLS_PATH} faidx {input} {wildcards.interval} -o {output}
"""
rule run_TRF:
    input:
        "temp/{interval}.fasta"
    output:
        dat = temp("temp/{interval}.dat")
    threads:
        4
    params:
        dat_path = (temp("temp/{interval}.dat"))
    shell:"""
fasta=$(realpath {input})
fasta1=$(realpath {params.dat_path})
cd temp/
pwd
echo {params.dat_path}
echo ${{fasta1}}
{TRF_PATH} ${{fasta}} 2 7 7 80 10 50 500 -h -l 20 -ngs > ${{fasta1}}
cd ..
"""
rule trf_to_bed:
	input:
		dats = expand(rules.run_TRF.output.dat, interval=INTERVALS)
	output:
		bed = f"{PREFIX}.trf.bed"
	threads:
		1
	run:
		trf = []
		header = '#chr start end PeriodSize CopyNumber ConsensusSize PercentMatches PercentIndels Score A C G T Entropy Motif Sequence'.split()
		for datf in input.dats:
			chrom = None
			sys.stderr.write( "\r" + datf )			
			with open(datf, 'r') as dat:
				for line in dat:
					splitline = line.split()
					if( line.startswith("Sequence:") ):
						chrom = int(line.split()[1].strip())
						#sys.stderr.write(chrom + "\n")
					elif( line.startswith("@") ):
						chrom = splitline[0][1:].strip() # grab everything after the @ in the first word
					else:
						# Catch index errors when line is blank
						try:
							# Check if in header sequence (all non-header lines start with an int: start pos)
							try:
								int(splitline[0])
							except ValueError:
								continue
							trf.append([chrom] + splitline[ 0: (len(header)-1) ] )
						except IndexError:
							pass
		trf = pd.DataFrame(trf, columns=header)
		print(trf.shape)
		
		trf["start"] = trf["start"].astype(int) - 1
		trf.sort_values(by=["#chr", "start"], inplace=True)
		print("done sorting trf")

		trf.to_csv(output.bed, sep="\t", index=False)
rule process:
	input:
		file = f"{PREFIX}.trf.bed"
	output:
		telomere = f"{PREFIX}telo.txt"
	shell:"""
	(grep -E 'CCCTAA' {input} | awk '{{
    count = 0  # 重置计数器，以便在处理新的行时重新计数
    for(j=1; j<=length($16); j++) {{
        if (substr($16, j, 6) == "CCCTAA") {{
            count++
        }}
    }}
    if (count>30){{print $1,$2,$3,count}}
}}') > {output}

(grep -E 'TTAGGG' {input} | awk '{{
    count = 0  # 重置计数器，以便在处理新的行时重新计数
    for(j=1; j<=length($16); j++) {{
        if (substr($16, j, 6) == "TTAGGG") {{
            count++
        }}
    }}
    if (count>30){{print $1,$2,$3,count}}
}}') >> {output}
	"""