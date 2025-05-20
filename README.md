# LSGvar -- Large-Scale Genomic VARiation caller
LSGvar is a caller for comprehensive large-scale structural variants detection based on assemblies.
## Configuration
Install LSGvar through conda:
```shell
##Before using LSGvar, configure the corresponding environment first.
git clone https://github.com/Hanjunmin/LSGvar.git
cd LSGvar-main
conda env create -f environment.yml 
```
*If you cannot use conda, you can obtain the environment dependencies by pulling the docker image:
```shell
docker pull crpi-2hir7evq6wnxr5no.cn-shanghai.personal.cr.aliyuncs.com/feifeizhou/lsgvar:v0.1
```

## Parameters
```shell
LSGVAR --help
 _        _____    _____  __       __   _       ______
| |      / ____|  / ____| \ \     / /  / \     |  __  \
| |     | (___   | |  __   \ \   / /  / - \    | |__| |
| |      \___ \  | | |_ \   \ \ / /  / /_\ \   |  __  /
| |___   ____) | | |___) |   \ - /  / /___\ \  | |  \ \
|_____| |_____/   \_____/     \_/  /_/     \_\ |_|   \_\

Large-Scale Genetic VARiation caller


Run command:
LSGVAR -r ref.fa -q1 hap1.fa -q2 hap2.fa -p1 hap1.paf -p2 hap2.paf -cp1 PTR_hap1_pair.tsv -cp2 PTR_hap2_pair.tsv -cen chm13_cen.tsv -telo chm13_telo.tsv -m cts

usage: LSGVAR [-h] -r REF -q1 HAP1 [-q2 HAP2] [-p1 PAF1] [-p2 PAF2] -cp1 PAIRS1 [-cp2 PAIRS2] [-c CLUSTER] [-d DELLENGTH] [-i INVCLUSTER] [-cen CENTROMERE] [-telo TELOMERE] -m
              {ctn,cts} [-vt VARIANT] [--nosnv]

options:
  -h, --help            show this help message and exit

Input Files:
  -r REF, --ref REF     Reference genome for variants calling
  -q1 HAP1, --hap1 HAP1
                        One query genome (Which is one haplotype of one species genome)
  -q2 HAP2, --hap2 HAP2
                        Another query genome (Which is another haplotype of the species genome)
  -p1 PAF1, --paf1 PAF1
                        Alignment of hap1 (Which contains the CIGAR infomation)
  -p2 PAF2, --paf2 PAF2
                        Alignment for hap2 (Which contains the CIGAR infomation)
  -cp1 PAIRS1, --pairs1 PAIRS1
                        Chromsome pairs of query genome (hap1) and reference
  -cp2 PAIRS2, --pairs2 PAIRS2
                        Chromsome pairs of query genome (hap2) and reference
  -c CLUSTER, --cluster CLUSTER
                        Clustering parameter for filtering, where a smaller value results in a stricter filter [200000]
  -d DELLENGTH, --dellength DELLENGTH
                        A desired deletion length for alignments, where a larger value enforces a stricter filter [300000]
  -i INVCLUSTER, --invcluster INVCLUSTER
                        Clustering parameter for inversion calling [700000]
  -cen CENTROMERE, --centromere CENTROMERE
                        A centromere file which is used to filter out the alignment of complex regions that may not be well aligned [False]
  -telo TELOMERE, --telomere TELOMERE
                        A telomere file which is used to filter out the alignment of complex regions that may not be well aligned [False]
  -m {ctn,cts}, --mode {ctn,cts}
                        Analysis mode: ctn (do not remove centromere and telomere alignments) or cts (remove) [ctn]

Additional arguments:
  -vt VARIANT, --variant VARIANT
                        Comma-separated variant types to generate final result (default: all). Options: snv, ins, del, inv, trans, sdr, dup, highdup
  --nosnv               Params for skipping SNV/Indel identification and merge
```

## Usage
The test data is in the ${Tool_dir}/examples/, and genome/ contains three zipped sample genomic data, namely the chromosome 21 of chimpanzees (two haplotypes) and human (chr21.ptr.hap1.fa.gz, chr21.ptr.hap2.fa.gz and chr21.chm13.fa.gz), data/ contains the centromere, telomere annotation data of T2T-CHM13, and one homologous chromosome pair file, the align/ folder contains two paf files for aligning the query genome to the reference genome.
```shell
##After configuring LSGvar, you can create a new working folder in any path
mkdir LSGvar_work && cd LSGvar_work
##Activate the environment
conda activate LSGvar
##If you want to test this tool, change to the test genome directory and unzip them
cd ${Tool_dir}/examples/genome
gunzip chr21.chm13.fa.gz && gunzip chr21.ptr.hap1.fa.gz && gunzip chr21.ptr.hap2.fa.gz
##Then use the command to start SV identification
cd LSGvar_work

##If you don't have the alignment paf file, just provide the fasta path
LSGVAR -r ${Tool_dir}/examples/genome/chr21.chm13.fa -q1 ${Tool_dir}/examples/genome/chr21.ptr.hap1.fa -q2 ${Tool_dir}/examples/genome/chr21.ptr.hap2.fa -cp1 ${Tool_dir}/examples/data/PTR_hap1_pairs.tsv -cp2 ${Tool_dir}/examples/data/PTR_hap2_pairs.tsv -cen ${Tool_dir}/examples/data/chm13_cen.tsv -telo ${Tool_dir}/examples/data/chm13_telo.tsv -m cts

##If you already have the paf file, use the -p1 (-p2) parameters
LSGVAR -r ${Tool_dir}/examples/genome/chr21.chm13.fa -q1 ${Tool_dir}/examples/genome/chr21.ptr.hap1.fa -q2 ${Tool_dir}/examples/genome/chr21.ptr.hap2.fa -p1 ${Tool_dir}/examples/align/align_hap1.paf -p2 ${Tool_dir}/examples/align/align_hap2.paf -cp1 ${Tool_dir}/examples/data/PTR_hap1_pairs.tsv -cp2 ${Tool_dir}/examples/data/PTR_hap2_pairs.tsv -cen ${Tool_dir}/examples/data/chm13_cen.tsv -telo ${Tool_dir}/examples/data/chm13_telo.tsv -m cts
```

## SV-annotation：

The final result can be found in `${work_dir}/results`.

**LSGvarhap1(2).bed**:

|Label     |annotation                                                |
| ----------------- | ------------------------------------------------------------ |
| SNV  | SNV |
| DEL          | Deletion (<50bp, >=50bp)                                              |
| INS           | Insertion (<50bp, >=50bp)                                                |
| DUP           | Duplication                                              |
| INV           | Inversion                                              |
| SDR           | Structure Divergent Reigions                                         |
| TRANS           | Translocation                                              |

**hap1(hap2)cigarsdr.vcf**:

INS、DEL、SNV、INV

**sortLSGvarall.vcf.gz**(hap1+hap2)

**LSGvar.bed**(hap1+hap2)
