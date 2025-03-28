import os
import sys
import click
import math
from pathlib import Path
from pie.module.sort_chrom import sort_chrom
import logging

@click.command()
@click.option("-i", "--input", required = True, type = str, help = "Input vcf/vcf.gz file for evaluation")
@click.option("-n", "--name", default = "Sample", type = str, help = "User defined sample name, [default: Sample]")
@click.option("-c", "--compare", required = True, type = str, help = "Ground truth vcf/vcf.gz file for comparison")
@click.option("-r", "--ref", required = True, type = str, help = "Reference file fasta file (.fasta or .fa) or fasta index file (.fai)")
@click.option("-o", "--output", required = True, type = str, help = "Output file prefix,  such as -o ./test/output_name")
@click.option("-t", "--threads", default = 24, type = int, help = "Maximum numbers of parallel threads [default: 24]")
@click.option("-m", "--max-len", default = 250 * 10**3, type = int, help = "Maximum variant distance for pairwise calculation [default: 250000]")
@click.option("-b", "--bed", default = None, type = str, help = r".bed file specifying genomic regions to include [default: None]")
@click.option("--min-sv", default = 30, type = int, help = "Minimal length threshold of Structral Variant [default: 30, ALT length > 30 bp is SV]")
@click.option("--chrom", default = ','.join(["chr" + str(i) for i in range(1, 23)]), type = str, help = "Chromosome to evaluate,use comma to join chromosome name e.g. --chrom chr1,chr2,chr3 [default:chr1,chr2,chr3,...,chr22]")
@click.option("--sexchrom", default = "chrX,chrY", type = str, help = "Sex chromosme,use comma to join chromosome name e.g. --sexchrom chrX,chrY [default: chrX,chrY]")
@click.option("--mincount", default = 2, type = int, help = "Minimum numbers of phased sites in a phase block [default: 2]")
@click.option("--block", is_flag = True, help = r"Output phasing block start and end positions in a BED file")
@click.option("--no-sex", is_flag = True, help = "Ignore sex chromosome")
@click.option("--canonical", is_flag = True, help = "Canonical mode, only evaluate single mutation SNV ignore double heterozygous site")
@click.option("--only-snv", is_flag = True, help = "Only evaluate single nucleotide variation")
@click.option("--only-indel", is_flag = True, help = "Only evaluate insertion and deletion")
@click.option("--only-sv", is_flag = True, help = "Only evaluate structural variant")
@click.option("--no-snv", is_flag = True, help = "Ignore single nucleotide variation")
@click.option("--no-indel", is_flag = True, help = "Ignore insertion and deletion")
@click.option("--no-sv", is_flag = True, help = "Ignore structural variant")
@click.option("--no-double", is_flag = True, help = "Ignore double heterozygous site")
@click.option("--no-sort", is_flag = True, help = "Do not sort chromosome or regions, directly use the input order")
@click.option("--verbose", is_flag = True, help = "Enable verbose mode, printing parameters and progress to standard output")
@click.version_option(version="0.7.0", prog_name = r"Phasing all-in-one evaluator (pie), based on Python 3.7+")
def cli(input, name, compare, ref, output, threads, max_len, bed, min_sv, chrom, sexchrom, mincount, canonical, block, no_sex, only_snv, only_indel, only_sv, no_snv, no_indel, no_sv, no_double, no_sort, verbose) -> tuple:

    # set logging
    logging.basicConfig(level = logging.DEBUG, format = "%(asctime)s - %(levelname)s - %(message)s")

    # clean parameters
    input = os.path.abspath(input)
    compare = os.path.abspath(compare)
    ref = os.path.abspath(ref)
    
    if bed:
        bed = os.path.abspath(bed)
    
    if canonical:
        only_snv = True
        no_double = True
    

    # check parameters
    input_vcf = Path(input)
    compare_vcf = Path(compare)
    if not (input_vcf.suffix == ".vcf" or input_vcf.suffix == ".bcf" or input_vcf.suffixes == [".vcf", ".gz"]):
        raise ValueError(r"-i or --input must be .vcf or .vcf.gz file. You provided a file with a different extension")
    if not (compare_vcf.suffix == ".vcf" or compare_vcf.suffix == ".bcf" or compare_vcf.suffixes == [".vcf", ".gz"]):
        raise ValueError(r"-c or --compare must be .vcf or .vcf.gz file. You provided a file with a different extension")     

    ref_file = Path(ref)
    if not (ref_file.suffix == ".fa" or ref_file.suffix == ".fasta" or ref_file.suffix == ".fai"):
        raise ValueError(r"-r or --ref must be .fa, .fasta, or .fai file. You provided a file with a different extension") 

    if bed:
        bed_file = Path(bed)
        if not (bed_file.suffix == ".bed"):
            raise ValueError(r"--bed must be .bed file. You provided a file with a different extension")

    output_dir = os.path.dirname(output)
    if not (os.path.exists(output_dir)):
        raise ValueError(f"The output directory {output_dir} does not exist")
    else:
        if not os.path.isdir(output_dir):
            raise ValueError(f"The output directory {output_dir} exists but is not a directory")        

    if threads > os.cpu_count():
        logging.warning(f"-t or --threads set threads exceed the system CPU thread count ({os.cpu_count()})")
        

    # process chromosome area
    chrom = [i.strip() for i in chrom.split(',')]
    if not chrom:
        raise ValueError("-c or --chrom must be in comma joined string e.g. chr1,chr2,ch3")
    
    sexchrom = [i.strip() for i in sexchrom.split(',')]
    if not chrom:
        raise ValueError("--sexchrom must be in comma joined string e.g. chr1,chr2,ch3")

    if no_sex:
        clean_chrom = []
        for i in chrom:
            if i not in sexchrom:
                clean_chrom.append(i)
        chrom = clean_chrom
    
    # sort chromosome
    if not no_sort:
        chrom = sort_chrom(chrom)

    # print parameters in verbose mode
    if verbose:
        ctx = click.get_current_context()
        logging.info("Command parameters:")
        for param in ctx.command.params:
            param_name = param.name
            if param_name != "version":
                param_value = ctx.params[param_name]
                logging.info(f"{param_name}: {param_value}")

    return (input, name, compare, ref, output, threads, max_len, bed, min_sv, chrom, sexchrom, mincount, canonical, block, no_sex, only_snv, only_indel, only_sv, no_snv, no_indel, no_sv, no_double, no_sort, verbose)


