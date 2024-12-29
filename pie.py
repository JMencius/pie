import click
import sys
import os
import time
import logging
from multiprocessing import Pool
from scripts.parse_vcf import get_sample_name
from scripts.intersect import intersect
from scripts.evaluation import blockwise_evaluate
from scripts.cal_ref import get_ref_len
from scripts.overall_metrics import cal_NGx0
from scripts.cal_truth_metrics import cal_truth_pairs
from scripts.F1_related import cal_F1_related
from multiprocessing import Pool
from scripts.sort_key import sort_key
from IO.read_bed import read_bed
from IO.read_vcf import read_vcf
from IO.write_results import write_results
from IO.write_stats import write_stats
from process.prefilter import prefilter


PWD = os.path.dirname(os.path.realpath(__file__))

@click.command()
@click.option("-i", "--input", required = True, help = "Input vcf/bcf file for evaluation")
@click.option("-n", "--name", default = None, help = "User defined sample name of the input vcf file, [default: `extract from vcf`]")
@click.option("-c", "--compare", required = True, help = "Truth vcf/bcf file for comparsion")
@click.option("-r", "--ref", required = True, help = "Reference file fasta or fasta.fai")
@click.option("-o", "--output", required = True, help = "Output tsv file prefix, path can be added before the prefix, such as -o /test/output_name")
@click.option("-t", "--threads", default = 24, help = "Maximum numbers of parallel threads")
@click.option("--bed", default = None, help = r"Regions to only include, defined in bed file")
@click.option("--min-sv", default = 30, help = "Minimal length of Structral Variant")
@click.option("--chrom", default = ','.join([str(i) for i in range(1, 23)] + ['X', 'Y']), help = "Chromosome to evaluate,use comma to join e.g. --chrom 1,2,3 [default:1-23, X, Y]")
@click.option("--mincount", default = 2, help = "Minimum numbers of phased sites in a phase block [default: 2]")
@click.option("--no-sex", is_flag = True, default = True, help = "Ignore sex chromosome")
@click.option("--canonical", is_flag = True, help = "Canonical mode, only evaluate single mutation SNV")
@click.option("--only-snv", is_flag = True, help = "Only evaluate Single Nucleotide Variation")
@click.option("--no-sv", is_flag = True, help = "Ignore Structural Variant")
@click.option("--no-indel", is_flag = True, help = "Ignore insertion and deletion")
@click.option("--no-double", is_flag = True, help = "Ignore double heterozygous site")
@click.option("--verbose", is_flag = True, help = "Verbose mode print intermediate results to stdout")
@click.version_option(version="es-0.4.0", prog_name = r"phasing all-in-one evaluator(pie), based on Python 3.7+")
def main(input, name, compare, ref, output, threads, bed, min_sv, chrom, mincount, canonical, no_sex, only_snv, no_sv, no_indel, no_double, verbose):
    if verbose:
        start_time = time.time()
    
    logging.basicConfig(level = logging.DEBUG, format = "%(asctime)s - %(levelname)s - %(message)s")

    # get aboslute path
    logging.info(f"Processing parameters")
    input = os.path.abspath(input)
    compare = os.path.abspath(compare)
    ref = os.path.abspath(ref)
        
    if bed:
        bed = os.path.abspath(bed)
    
    if canonical:
        only_snv = True
        no_double = True
    
    # process chromosome area
    chrom = [i.upper() for i in chrom.split(',')]
    if no_sex:
        if 'X' in chrom:
            chrom.remove('X')
        if 'Y' in chrom:
            chrom.remove('Y')
    chrom.sort(key = sort_key)

    # print parameters in verbose mode
    if verbose:
        ctx = click.get_current_context()
        logging.info("Command parameters:")
        for param in ctx.command.params:
            param_name = param.name
            if param_name != "version":
                if param_name == "chrom":
                    logging.info(f"chrom: {chrom}")
                    continue
                param_value = ctx.params[param_name]
                logging.info(f"{param_name}: {param_value}")
    

    # read reference .fa or .fai
    logging.info("Parsing reference file")
    len_dict = get_ref_len(ref, chrom)
    total_ref_len = sum(len_dict.values())


    # get target area
    if bed:
        target = read_bed(bed)
    
    
    # read and filter two vcf files
    logging.info("Reading query VCF file")
    with Pool(threads) as p:
        query_vcf = p.starmap(read_vcf, [(input, c) for c in chrom])

    logging.info("Reading truth VCF file")
    with Pool(threads) as p:
        truth_vcf = p.starmap(read_vcf, [(compare, c) for c in chrom])
    ##print(list(truth_vcf[0].items())[:10])

    
    logging.info("Evaluate genotype and extract phased")
    with Pool(threads) as q:
        prefiltered = q.starmap(prefilter, [(query_vcf[t], truth_vcf[t], chrom[t], min_sv, only_snv, no_sv, no_indel, no_double) for t in range(len(chrom))])
    logging.info("Finish evaluate genotype and extract phased")


    logging.info("Writing variant stat file")
    write_stats([i[2] for i in prefiltered], output ,chrom)


    # calculate truth metrics
    logging.info("Calculating truth metrics")

    with Pool(threads) as t:
        pairs_results = t.map(cal_truth_pairs, [prefiltered[i][1] for i in range(len(chrom))])


    logging.info("Intersecting two vcf files")
    with Pool(threads) as q:
        intersect_results = q.starmap(intersect, [(prefiltered[t][0], prefiltered[t][1], chrom[t], mincount, min_sv) for t in range(len(chrom))])
        #intersect_results = q.starmap(intersect, [(query_vcf[t_chr], truth_vcf[t_chr], t_chr, mincount, min_sv) for t_chr in chrom])
    
    logging.info("Evluating blocks")
    with Pool(threads) as r:
        evaluation_results = r.starmap(blockwise_evaluate, [(intersect_results[i], i, len_dict, verbose) for i in range(len(intersect_results))])
    
    logging.info("Calculating NG50 and NG90")
    phase_len = list()
    for eva in evaluation_results:
        phase_len += eva["length_list"]
    phase_len.sort(reverse = True)
    NG50 = cal_NGx0(phase_len, total_ref_len, 50)
    NG90 = cal_NGx0(phase_len, total_ref_len, 90)
    
    # calculate F1-score
    evaluation_results = cal_F1_related(pairs_results, evaluation_results)
    
    # Output to files
    if (not name):
        name = get_sample_name(input)
    write_results(evaluation_results, chrom, (NG50, NG90), output, name, verbose)
    logging.info("ALL DONE")

    if verbose:
        end_time = time.time()
        logging.info(f"Total processing time is {end_time - start_time} seconds.")



if __name__ == "__main__":
    main()

