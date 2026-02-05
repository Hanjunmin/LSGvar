# MiTool
# LSGvar -- Large-Scale Genomic VARiation caller
LSGvar is a caller for comprehensive large-scale structural variants detection based on assemblies, demonstrating superior performance in cross-species variant identification, particularly in inversion detection, compared to existing tools.
## Configuration
Install LSGvar through conda:
```shell
##Before using LSGvar, configure the corresponding environment first. You can build a new env before running.
git clone https://github.com/Hanjunmin/LSGvar.git
cd LSGvar-main
conda env create -f environment.yml 
```
*If you cannot use conda, you can obtain the environment dependencies by pulling the docker or singularity image:
```shell
docker pull crpi-2hir7evq6wnxr5no.cn-shanghai.personal.cr.aliyuncs.com/feifeizhou/lsgvar:v0.1
##then enter this container like this:
docker run -it -v $(pwd):/LSGVAR lsgvar:v0.1

singularity pull --arch amd64 library://feifeizhou/tool/lsgvar:v0.1
##The tool is at /LSGVAR/main/LSGVAR path under the container, then you can use it like this:
singularity exec lsgvar_v0.1.sif /LSGVAR/main/LSGVAR -r ref.fa -q1 query1.fa -cp1 pair1.tsv -m ctn
```

## Parameters
```shell
LSGVAR --help
╔════════════════════════════════════════════════════════╗
║                                                        ║
║   ██╗     ███████╗ ██████╗ ██╗   ██╗ █████╗ ██████╗    ║
║   ██║     ██╔════╝██╔════╝ ██║   ██║██╔══██╗██╔══██╗   ║
║   ██║     ███████╗██║  ███╗██║   ██║███████║██████╔╝   ║
║   ██║     ╚════██║██║   ██║██║   ██║██╔══██║██╔══██╗   ║
║   ███████╗███████║╚██████╔╝╚██████╔╝██║  ██║██║  ██║   ║
║   ╚══════╝╚══════╝ ╚═════╝  ╚═════╝ ╚═╝  ╚═╝╚═╝  ╚═╝   ║
║                                                        ║
║                       L S G V A R                      ║
╚════════════════════════════════════════════════════════╝
Large-Scale Genetic VARiation caller


Run command:
LSGVAR -r ref.fa -q1 hap1.fa -q2 hap2.fa -p1 hap1.paf -p2 hap2.paf -cp1 PTR_hap1_pair.tsv -cp2 PTR_hap2_pair.tsv -cen chm13_cen.tsv -telo chm13_telo.tsv -m cts -s PTR

usage: test_LSGVAR [-h] -r REF -q1 HAP1 [-q2 HAP2] [-p1 PAF1] [-p2 PAF2] -cp1 PAIRS1 [-cp2 PAIRS2] [-c CLUSTER] [-dl DELLENGTH]
                   [-inv INVCLUSTER] [-cen CENTROMERE] [-telo TELOMERE] -m {ctn,cts} [-d DISTANCE] [-t THREADS] [-k CHUNK_SIZE]
                   [-s SAMPLE_NAME] [-mdist MAX_DISTANCE] [-sdist SMALL_DISTANCE] [-sim SIMILARITY_THRESHOLD]

options:
  -h, --help            show this help message and exit

Input Files:
  -r REF, --ref REF     Reference genome for variants calling
  -q1 HAP1, --hap1 HAP1
                        One query genome (Which is one haplotype of one species genome and needs to be scaffolded to chromosome
                        level using RagTag to ensure better genome quality for more reliable variant calling results.)
  -q2 HAP2, --hap2 HAP2
                        Another query genome (Which is another haplotype of the species genome and needs to be scaffolded to
                        chromosome level using RagTag to ensure better genome quality for more reliable variant calling results.)
  -p1 PAF1, --paf1 PAF1
                        Alignment of haplotype1 (Which contains the CIGAR infomation) [Recommend mapping tool: minimap2].
  -p2 PAF2, --paf2 PAF2
                        Alignment for haplotype2 (Which contains the CIGAR infomation) [Recommend mapping tool: minimap2].
  -cp1 PAIRS1, --pairs1 PAIRS1
                        Homologous chromsome pairs of query genome (hap1) and reference.
  -cp2 PAIRS2, --pairs2 PAIRS2
                        Homologous chromsome pairs of query genome (hap2) and reference.
  -c CLUSTER, --cluster CLUSTER
                        Clustering parameter for filtering chaos alignments, where a smaller value results in a stricter filter
                        [200000].
  -dl DELLENGTH, --dellength DELLENGTH
                        A desired deletion length for alignments, where a larger value enforces a stricter filter [300000].
  -inv INVCLUSTER, --invcluster INVCLUSTER
                        Clustering parameter for inversion calling [700000].
  -cen CENTROMERE, --centromere CENTROMERE
                        A centromere file which is used to filter out the alignment of complex regions that may not be well aligned
                        [False].
  -telo TELOMERE, --telomere TELOMERE
                        A telomere file which is used to filter out the alignment of complex regions that may not be well aligned
                        [False].
  -m {ctn,cts}, --mode {ctn,cts}
                        Analysis mode: ctn (do not remove centromere and telomere alignments) or cts (remove) [ctn].
  -d DISTANCE, --distance DISTANCE
                        Parameters used to identify INS and DEL: Variations where the distance between the start and end positions
                        of the REF or QUERY is less than (or equal with) d will be identified as SVs. A lower value indicates
                        stricter criteria for SV identification [0bp].
  -t THREADS, --threads THREADS
                        Multi threads for inversion re-identification from SDR.
  -k CHUNK_SIZE, --chunk_size CHUNK_SIZE
                        Process line of each thread in inversion re-identification.
  -s SAMPLE_NAME, --sample_name SAMPLE_NAME
                        Sample name used to generate vcf.

Additional arguments:
  -mdist MAX_DISTANCE, --max_distance MAX_DISTANCE
                        Max reference distance for two allele to merge [500bp].
  -sdist SMALL_DISTANCE, --small_distance SMALL_DISTANCE
                        Max reference distance for SSV (small variants) merge [10bp].
  -sim SIMILARITY_THRESHOLD, --similarity_threshold SIMILARITY_THRESHOLD
                        The similarity of variants, used to merge the variants of two haplotypes [0.8].
```

## Usage  
### About query genome
The scaffold assembly may contain chromosomes that have over 50% aligned in reverse orientation compared to the reference genome. The script ``LSGvar/scripts/rc_chrom.sh`` automatically detects such cases by analyzing alignment data (in paf format) generated from Minimap2. We recommend that you first align the query genome to the reference genome to check for this situation, and then perform reverse complementation to identify more accurate inversions.

### Test
The test data is in ``${Tool_dir}/examples/``  

The ``genome/`` folder contains three zipped sample genomic data, namely the chromosome 21 of chimpanzees (two haplotypes) and human:  
- ``chr21.ptr.hap1.fa.gz``  
- ``chr21.ptr.hap2.fa.gz``  
- ``chr21.chm13.fa.gz``  

The ``data/`` folder contains:  
- Centromere and telomere annotation data of T2T-CHM13  
- Homologous chromosome pair files  

The ``align/`` folder contains two PAF files for aligning the query genome to the reference genome.  

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

## Results：

The final result can be found in `${work_dir}/results`.

**sortLSGvar_all.vcf**(INDEL、SNV、SV of two haplotypes)

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
| INV-INV | Nested inversion |

**hap1(hap2)cigarsdr.vcf**:

INS、DEL、SNV、INV

**LSGvar.bed**(hap1+hap2)

## Getting help
If you have any questions about its use, please raise an issue here.

## Citing LSGVAR
If you use LSGvar in your research, please cite:
https://github.com/Hanjunmin/LSGvar/
> To be updated
