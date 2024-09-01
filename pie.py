import click
import sys
import os
import time
from scripts.parse_vcf import process_filter_vcf
from multiprocessing import Pool



@click.command()
@click.option("-i", "--input", required = True, help = "Input vcf/bcf file for evaluation")
@click.option("-c", "--compare", required = True, help = "Truth vcf/bcf file for comparsion")
@click.option("-r", "--ref", default = "GCA_000001405.15_GRCh38_no_alt", help = "Reference file [default: GCA_000001405.15_GRCh38_no_alt]")
@click.option("-o", "--output", default = "Stdout", help = r"Output [default: Stdout]")
@click.option("-t", "--threads", default = 12, help = "Maximum numbers of parallel threads")
@click.option("--fbed", default = None, help = r"Bed file to filter out, such as centromere region in data/hg38_centromere.bed")
@click.option("--min-sv", default = 30, help = "Minimal length of Structral Variant")
@click.option("--chrom", default = ','.join([str(i) for i in range(1, 23)] + ['X', 'Y']), help = "Chromosome to evaluate,use comma to connect e.g. --chrom 1,2,3")
@click.option("--mincount", default = 2, help = "Minimum numbers of phased sites in a phase block")
@click.option("--no-sex", is_flag = True, help = "Ignore sex chromosome")
@click.option("--canonical", is_flag = True, help = "Canonical mode, only evaluate single mutation SNV")
@click.option("--only-snv", is_flag = True, help = "Only evaluate Single Nucleotide Variation")
@click.option("--no-sv", is_flag = True, help = "Ignore Structural Variant")
@click.option("--no-indel", is_flag = True, help = "Ignore insertion and deletion")
@click.option("--no-double", is_flag = True, help = "Ignore double heterozygous site")
@click.option("--no-centro", is_flag = True, help = "Ignore centromere region")
@click.option("--only-num", is_flag = True, help = "Export only numbers")
@click.option("--verbose", is_flag = True, help = "Verbose mode print some results to stdout")
@click.version_option(version="0.1.0", prog_name = r"phasing all-in-one evaluator(pie), based on Python 3.7+")
def main(input, compare, ref, output, threads, fbed, min_sv, chrom, mincount, canonical, no_sex, only_snv, no_sv, no_indel, no_double, no_centro, only_num, verbose):
    start_time = time.time()
    
    # get aboslute path
    input = os.path.abspath(input)
    compare = os.path.abspath(compare)
    if ref:
        ref = os.path.abspath(ref)
    output = os.path.abspath(output)
    
    # process chromosome area
    chrom = set(chrom.split(','))

    # print parameters
    ctx = click.get_current_context()
    click.echo('Command parameters:')
    
    for param in ctx.command.params:
        param_name = param.name
        if param_name != "version":
            param_value = ctx.params[param_name]
            click.echo(f'{param_name}: {param_value}')
    print('\n')
    
    ###sys.exit(0)

    # simultaneously read two vcf files
    print("Start reading vcf files")
    read_threads = min(2, threads)
    with Pool(read_threads) as p:
    ##process_filter_vcf(filename: str, fbed: str, min_sv: int, chrom: set, no_sex : bool, canonical : bool, only_snv : bool, no_sv : bool, no_indel : bool, no_double : bool, no_centro : bool)
        query_vcf, truth_vcf = p.starmap(process_filter_vcf, [(v, fbed, min_sv, chrom, no_sex, canonical, only_snv, no_sv, no_indel, no_double, no_centro) for v in [input, compare]])
    
    print("Query vcf phase count: ", query_vcf[1])
    print("Query vcf unphase count: ", query_vcf[2])
    print("Truth vcf phase count: ", truth_vcf[1])
    print("Truth vcf unphase count: ", truth_vcf[2])

    end_time = time.time()
    print(f"ALL DONE")
    print(f"Total processing time is {end_time - start_time} seconds.")



if __name__ == "__main__":
    main()
