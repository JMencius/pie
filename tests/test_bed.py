import subprocess
import os

script_dir = os.path.dirname(os.path.abspath(__file__))


def test_bed_specify_region():
    result = subprocess.run(
        ["pie", "-i", f"{script_dir}/small_query.vcf", 
        "-c", f"{script_dir}/small_truth.vcf",
        "-r", f"{script_dir}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.fai",
        "--bed", f"{script_dir}/chr1_region.bed",
        "-o", f"{script_dir}/chr1_region"],
        capture_output = True,
        text = True             
    )
    
    assert result.returncode == 0, f"Command failed with error: {result.stderr}"


