import click
import sys
import os
import time
import tracemalloc
from scripts.parse_vcf import process_filter_vcf
from scripts.parse_vcf import get_sample_name
from scripts.intersect import intersect
from scripts.evaluation import blockwise_evaluate
from scripts.write_results import write_results
from scripts.cal_ref import cal_total_ref
from scripts.overall_metrics import cal_NGx0
from multiprocessing import Pool
from scripts.sort_key import sort_key

PWD = os.path.dirname(os.path.realpath(__file__))

@click.command()
@click.option("-i", "--input", required = True, help = "Input vcf/bcf file for evaluation")
@click.option("-n", "--name", default = None, help = "User defined sample name of the input vcf file, [default: `extract from vcf`]")
@click.option("-c", "--compare", required = True, help = "Truth vcf/bcf file for comparsion")
@click.option("-r", "--ref", default = f"{os.path.join(PWD, 'ref', 'GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta')}", help = "Reference file [default:GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta]")
@click.option("-o", "--output", required = True, help = "Output tsv file prefix, path can be added before the prefix, such as -o /test/output_name")
@click.option("-t", "--threads", default = 24, help = "Maximum numbers of parallel threads")
@click.option("--fbed", default = None, help = r"Bed file to filter out, such as centromere region in data/hg38_centromere.bed")
@click.option("--min-sv", default = 30, help = "Minimal length of Structral Variant")
@click.option("--chrom", default = ','.join([str(i) for i in range(1, 23)] + ['X', 'Y']), help = "Chromosome to evaluate,use comma to connect e.g. --chrom 1,2,3 [default:1-23, X, Y]")
@click.option("--mincount", default = 2, help = "Minimum numbers of phased sites in a phase block [default: 2]")
@click.option("--no-sex", is_flag = True, help = "Ignore sex chromosome")
@click.option("--canonical", is_flag = True, help = "Canonical mode, only evaluate single mutation SNV")
@click.option("--only-snv", is_flag = True, help = "Only evaluate Single Nucleotide Variation")
@click.option("--no-sv", is_flag = True, help = "Ignore Structural Variant")
@click.option("--no-indel", is_flag = True, help = "Ignore insertion and deletion")
@click.option("--no-double", is_flag = True, help = "Ignore double heterozygous site")
@click.option("--no-centro", is_flag = True, help = "Ignore centromere region")
@click.option("--verbose", is_flag = True, help = "Verbose mode print intermediate results to stdout")
@click.option("--test", is_flag = True, help = "Run test sample")
@click.version_option(version="es-0.1.1", prog_name = r"phasing all-in-one evaluator(pie), based on Python 3.7+")
def main(input, name, compare, ref, output, threads, fbed, min_sv, chrom, mincount, canonical, no_sex, only_snv, no_sv, no_indel, no_double, no_centro, verbose, test):
    if verbose:
        start_time = time.time()
    
    # get aboslute path
    input = os.path.abspath(input)
    compare = os.path.abspath(compare)
    ref = os.path.abspath(ref)
    
    # process chromosome area
    chrom = [i.upper() for i in chrom.split(',')]
    if no_sex:
        if 'X' in chrom:
            chrom.remove('X')
        if 'Y' in chrom:
            chrom.remove('Y')
    chrom.sort(key = sort_key)
    if verbose:
        print(f"Current working chromosome is {chrom}")


    # print parameters in verbose mode
    if verbose:
        ctx = click.get_current_context()
        click.echo('Command parameters:')
        for param in ctx.command.params:
            param_name = param.name
            if param_name != "version":
                param_value = ctx.params[param_name]
                click.echo(f'{param_name}: {param_value}')


    # simultaneously read and filter two vcf files
    print("Start reading vcf files")
    read_threads = min(2, threads)
    with Pool(read_threads) as p:
        query_vcf, truth_vcf = p.starmap(process_filter_vcf, [(v, fbed, min_sv, chrom, no_sex, canonical, only_snv, no_sv, no_indel, no_double, no_centro) for v in [input, compare]])
    
    if verbose:
        print(r"Output format in list is [#single_count, #double_count]")
        print("Query vcf phase count: ", query_vcf[1])
        print("Query vcf unphase count: ", query_vcf[2])
        print("Truth vcf phase count: ", truth_vcf[1])
        print("Truth vcf unphase count: ", truth_vcf[2])
    
    print("Intersecting two vcf files")
    with Pool(threads) as q:
        intersect_results = q.starmap(intersect, [(query_vcf[0][t_chr], truth_vcf[0][t_chr], t_chr, mincount, min_sv) for t_chr in chrom])
    
    print("Evluating blocks")
    with Pool(threads) as r:
        evaluation_results = r.starmap(blockwise_evaluate, [(intersect_results[i], i, verbose) for i in range(len(intersect_results))])
    
    print("Calculating NG50 and NG90")
    total_ref_len = cal_total_ref(ref, chrom)
    phase_len = list()
    for eva in evaluation_results:
        phase_len += eva["length_list"]
    phase_len.sort(reverse = True)
    NG50 = cal_NGx0(phase_len, total_ref_len, 50)
    NG90 = cal_NGx0(phase_len, total_ref_len, 90)
    
    # Output to files
    if (not name):
        name = get_sample_name(input)
    write_results(evaluation_results, chrom, (NG50, NG90), output, name, verbose)
    print("ALL DONE")

    if verbose:
        end_time = time.time()
        print(f"Total processing time is {end_time - start_time} seconds.")



if __name__ == "__main__":
    main()

