import subprocess
import os

script_dir = os.path.dirname(os.path.abspath(__file__))


def test_standard_vcf_input():
    result = subprocess.run(
        ["pie", "-i", f"{script_dir}/small_query.vcf", 
        "-c", f"{script_dir}/small_truth.vcf",
        "-r", f"{script_dir}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.fai",
        "-o", f"{script_dir}/vcf_in"],
        capture_output = True,
        text = True             
    )
    
    assert result.returncode == 0, f"Command failed with error: {result.stderr}"

def test_vcf_gz_input():
    result = subprocess.run(
        ["pie", "-i", f"{script_dir}/small_query.vcf.gz", 
        "-c", f"{script_dir}/small_truth.vcf",
        "-r", f"{script_dir}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.fai",
        "-o", f"{script_dir}/vcf_gz_in"],
        capture_output = True,
        text = True             
    )
    
    assert result.returncode == 0, f"Command failed with error: {result.stderr}"


def test_bcf_input():
    result = subprocess.run(
        ["pie", "-i", f"{script_dir}/small_query.bcf", 
        "-c", f"{script_dir}/small_truth.vcf",
        "-r", f"{script_dir}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.fai",
        "-o", f"{script_dir}/bcf_in"],
        capture_output = True,
        text = True             
    )
    
    assert result.returncode == 0, f"Command failed with error: {result.stderr}"

def test_all_consistent():
    with open("vcf_in.overall.csv", 'r') as f1, open("vcf_gz_in.overall.csv", 'r') as f2, open("bcf_in.overall.csv", 'r') as f3:
        vcf_in = f1.readlines()
        vcf_gz_in = f2.readlines()
        bcf_in = f3.readlines()
    
        assert vcf_in == vcf_gz_in == bcf_in, "The results of three types of format input is not consistent"

