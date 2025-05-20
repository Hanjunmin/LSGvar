import os
import time
import subprocess
from pathlib import Path
import glob

def run_minimap(ref_path, hap_path, hap_label, chrom_pair_file):
    """1. Whole genome alignments"""
    cmd = [
        "minimap2",
        "-t", "12",
        "-cx", "asm20",
        "--secondary=no",
        "--eqx",
        "-K", "8G",
        "-s", "1000",
        ref_path,
        hap_path,
        "-o", f"align_{hap_label}.paf"
    ]
    
    print(f"1.Running minimap2 for {hap_label}...")
    
    result = subprocess.run(cmd, check=True)
    
    if result.returncode == 0:
        print(f"Minimap2 completed for {hap_label} at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}", flush=True)
    
    print(f"2.Filtering PAF for {hap_label}...")
    process_paf(f"align_{hap_label}.paf",chrom_pair_file,hap_label)

def process_paf(paf_file, chrom_pair_file, hap_label):
    """2. Filter extra alignments"""
    filter_script = Path(os.path.dirname(os.path.abspath(__file__))) / "LSGvar/pyFiles/one2multi_filter.py"
    p_c_chrlen_file = f"{hap_label}_chrlen.txt"
    output_file = f"align_{hap_label}.flt.paf"
    
    cmd = [
        "python", str(filter_script),
        "-m", chrom_pair_file,  ## create automately
        "-f", paf_file,
        "-1", "6",
        "-2", "1"
    ]
    
    print(f"2.Filtering PAF for {hap_label}...")
    
    with open(output_file, 'w') as f:
        result = subprocess.run(cmd, stdout=f, check=True)
        print(f"Filter finished at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}", flush=True)
    
    awk_cmd = [
        "awk", "{print $6,$7,$1,$2}",
        f"{output_file}"
    ]
    
    with open(p_c_chrlen_file, 'w') as f:
        subprocess.run(awk_cmd, stdout=f, check=True)

def process_filtering(hap_label, mode, cluster, dellength, centromere, telomere):
    """3. Filter centromere and telomere alignments"""
    nowdic = Path.cwd()
    
    saffire_dir = Path(f"saffire{hap_label}")
    saffire_dir.mkdir(exist_ok=True)
    p_c_chrlen_file = nowdic/f"{hap_label}_chrlen.txt"
    output_file = nowdic/f"{hap_label}_syntenic.tsv"
    
    # run scripts
    chaos_script = Path(os.path.dirname(os.path.abspath(__file__))) / "LSGvar/scripts/chaos_filt.r"
    
    if mode == "ctn":
        cmd = [
            "Rscript", str(chaos_script),
            str(p_c_chrlen_file), 
            f"{str(saffire_dir)}/", 
            f"align_{hap_label}.flt.paf", 
            f"align_{hap_label}.final.paf",
            str(cluster),
            str(dellength),
            str(mode)
        ]
        
        print(f"3.Running chaos filter for {hap_label} in ctn mode...")
        
        #result = subprocess.run(cmd, check=True)
        result = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        print(f"Chaos filter finished at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}", flush=True)
        
    elif mode == "cts":
        if centromere is None or telomere is None:
            print("Error: When using 'cts' mode, you must provide both --centromere and --telomere arguments.")
            sys.exit(1)
        cmd = [
            "Rscript", str(chaos_script),
            str(p_c_chrlen_file), 
            f"{str(saffire_dir)}/", 
            f"align_{hap_label}.flt.paf", 
            f"align_{hap_label}.final.paf",
            str(cluster),
            str(dellength),
            str(mode),
            str(centromere), str(telomere)
        ]
        
        print(f"3.Running chaos filter for {hap_label} in cts mode...")
        
        #result = subprocess.run(cmd, check=True)
        result = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        print(f"Chaos filter finished at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}", flush=True)
    
    cmd = f"awk 'BEGIN{{OFS=\"\\t\"}} {{print $6, $8, $9, ($8+$9)/2, $1, $3, $4, ($3+$4)/2, $5}}' align_{hap_label}.final.paf | sort -k1,1 -k2,2n > {output_file}"
    subprocess.run(cmd, shell=True, check=True)

def run_cluster_and_call(ref_path, hap_path, hap_label, invcluster):
    """4. Cluster and SV calling"""
    
    denSDR_dir = Path(f"denSDR{hap_label}")
    denSDR_dir.mkdir(exist_ok=True)
    sdrall_file = denSDR_dir/"SDRall.txt"
    
    cluster_script = Path(os.path.dirname(os.path.abspath(__file__))) / "LSGvar/scripts/denSDR.r"
    fun_script = Path(os.path.dirname(os.path.abspath(__file__))) / "LSGvar/scripts/denSDRfun.r"
    
    cmd = [
        "Rscript", str(cluster_script),
        str(fun_script),
        f"{hap_label}_syntenic.tsv",
        f"{hap_label}_chrlen.txt",
        f"{str(denSDR_dir)}/",
        str(invcluster)
    ]
    
    print(f"4.Running clustering for {hap_label}...")
    
    #result = subprocess.run(cmd, check=True)
    result = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    print(f"{hap_label} clustering finished at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}", flush=True)
    
    cmd = [
        "cat", f"{denSDR_dir}/*end.tsv", ">", str(sdrall_file)
    ]
    
    subprocess.run(" ".join(cmd), shell=True, check=True)
    
    # realign to find simple inversions from SDR
    realign_script = Path(os.path.dirname(os.path.abspath(__file__))) / "LSGvar/vcf/realign.sh"

    cmd = [
        "bash", str(realign_script),
        "-i", str(sdrall_file),
        "-r", ref_path,
        "-q", hap_path,
        "-o", str(denSDR_dir/"SDRall_final.txt")
    ]
    
    subprocess.run(cmd, check=True)
    print(f"{hap_label} realign finished at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}", flush=True)

def run_cigar_processing(ref_path, hap_path, hap_label):
    """ 5.extract variationas from CIGAR """
    cigar_script = Path(os.path.dirname(os.path.abspath(__file__))) / "LSGvar/pyFiles/CIGAR.py"
    temp_dir = Path("temp")
    temp_dir.mkdir(exist_ok=True)
    
    cmd = [
        "python", str(cigar_script),
        "--r", ref_path,
        "--q", hap_path,
        "--paf", f"align_{hap_label}.final.paf",
        "--o", f"{hap_label}cigar.txt"
    ]
    
    print(f"5.Generating CIGAR for {hap_label}...")
    
    subprocess.run(cmd, check=True)
    
    output_file = f"{hap_label}cigarend.txt"
    
    with open(output_file, 'w') as dest:
        with open(f"{hap_label}cigar.txt", 'r') as src:
            for line in src:
                dest.write(line)
        
        for file_path in glob.glob(os.path.join("temp", "*.cigar")):
            with open(file_path, 'r') as src:
                for line in src:
                    dest.write(line)

    cmd = [
        "find", str(temp_dir), "-type", "f", "-name", "*.cigar", "-delete"
    ]
    
    subprocess.run(cmd, check=True)
    print(f"CIGAR generated for {hap_label} finished at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}", flush=True)

def run_dup_filtering(hap_label):
    """6. dup filter"""
    dup_filt_script = Path(os.path.dirname(os.path.abspath(__file__))) / "LSGvar/vcf/dup_filt.sh"
    filt_script = Path(os.path.dirname(os.path.abspath(__file__))) / "LSGvar/pyFiles/dup_filt.py"
    
    paf_file = f"align_{hap_label}.final.paf"
    cigar_end_file = f"{hap_label}cigarend.txt"
    output_file = f"{hap_label}cigarout.txt"
    
    cmd = [
        "bash", str(dup_filt_script),
        paf_file,
        cigar_end_file,
        output_file,
        filt_script
    ]
    
    print(f"6.Running dup filter for {hap_label}...")
    
    subprocess.run(cmd, check=True)
    print(f"Dup filter for {hap_label} finished at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}", flush=True)


def generate_vcf(ref_path, hap_path, hap_label, variant_type):
    """7. Generate variation results -- vcf"""
    results_dir = Path(f"results")
    results_dir.mkdir(exist_ok=True)
    
    # 运行 cigar2vcf.sh
    cigar_to_vcf_script = Path(os.path.dirname(os.path.abspath(__file__))) / "LSGvar/vcf/cigar2vcf.sh"
    py_script = Path(os.path.dirname(os.path.abspath(__file__))) / "LSGvar/pyFiles/SDR_vcf.py"
    cigar_out_file = f"{hap_label}cigarout.txt"
    vcf_file = f"{results_dir}/{hap_label}cigarsdr.vcf"
    sdrall_final_file = f"denSDR{hap_label}/SDRall_final.txt"
    cigarsdr_txt_file = f"{hap_label}cigarsdr.txt"
    lsgvarend_bed_file = f"{results_dir}/LSGvar{hap_label}.bed"
 
    cmd = [
        "bash", str(cigar_to_vcf_script),
        cigar_out_file,
        vcf_file,
        sdrall_final_file,
        ref_path,
        hap_path,
        cigarsdr_txt_file,
        str(py_script),
        lsgvarend_bed_file,
        variant_type
    ]
    
    print(f"7.Generating VCF for {hap_label}...")
    subprocess.run(cmd, check=True)
    print(f"Generating VCF for {hap_label} finished at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}", flush=True)


def split_vcf(hap_label, variant_type):
    """8. split vcf"""
    split_script = Path(os.path.dirname(os.path.abspath(__file__))) / "LSGvar/vcf/splitfile.sh"
    hap_dir = Path(f"{hap_label}")
    hap_dir.mkdir(exist_ok=True)
    
    cmd = [
        "bash", str(split_script),
        f"{str(hap_dir)}/",
        f"results/{hap_label}cigarsdr.vcf",
        variant_type
    ]
    
    print(f"8.Splitting VCF for {hap_label}...")
    
    subprocess.run(cmd, check=True)
    print(f"Splitting VCF for {hap_label} finished at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}", flush=True)


def integrate_results(ref_path, hap1_dir, hap2_dir, variant_type):
    """9. merge two haplotypes' variation"""
    integrate_dir =  Path("results/integrate")
    integrate_dir.mkdir(exist_ok=True)

    cmd = f"less -S align_hap1.final.paf | awk 'BEGIN{{OFS=\"\\t\"}} {{print $6,$8,$9}}'| bedtools sort | bedtools merge -i - > results/integrate/hap1_paf.bed"
    subprocess.run(cmd, shell=True, check=True)
    
    cmd = f"less -S align_hap2.final.paf | awk 'BEGIN{{OFS=\"\\t\"}} {{print $6,$8,$9}}'| bedtools sort | bedtools merge -i - > results/integrate/hap2_paf.bed"
    subprocess.run(cmd, shell=True, check=True)

    cmd = [
        "bedtools", "intersect",
        "-a", f"{integrate_dir}/hap1_paf.bed",
        "-b", f"{integrate_dir}/hap2_paf.bed",
        "|", "bedtools", "sort",
        ">", f"{integrate_dir}/h1_h2intersec.bed"
    ]
    
    subprocess.run(" ".join(cmd), shell=True, check=True)
    
    phenotype_script = Path(os.path.dirname(os.path.abspath(__file__))) / "LSGvar/vcf/phenotype.sh"
    
    cmd = [
        "bash", str(phenotype_script),
        hap1_dir,
        hap2_dir,
        ref_path,
        variant_type
    ]
    
    subprocess.run(cmd, check=True)

    vcf2bedgt_script = Path(os.path.dirname(os.path.abspath(__file__))) / "LSGvar/vcf/vcf2bedGT.sh"
    
    cmd = [
        "bash", str(vcf2bedgt_script),
        f"results/sortLSGvarall.vcf.gz",
        f"results/LSGvarhap1.bed",
        f"results/LSGvarhap2.bed",
        f"results/LSGvar.bed",
        variant_type
    ]
    
    subprocess.run(cmd, check=True)
    print(f"LSGvar finished at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}")

def remove_temp_files():
    current_dir = Path(".")
    extensions = ('.txt', '.csv', '.fa', '.tsv')
    
    try:
        for file in current_dir.glob('*'):
            if file.suffix.lower() in extensions:
                file.unlink()  # Remove the file
    except Exception as e:
        sys.exit(1)
