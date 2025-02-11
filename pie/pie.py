import click
import sys
import os
import time
import logging
from multiprocessing import Pool
from module.sort_chrom import sort_chrom
from module.cal_ref import get_ref_len
from module.read_bed import read_bed
from module.read_vcf import read_vcf
from module.write_genotype import write_genotype
from module.evaluation import blockwise_evaluate


from scripts.parse_vcf import get_sample_name
from scripts.intersect import intersect
from scripts.evaluation import blockwise_evaluate
from scripts.overall_metrics import cal_NGx0
from scripts.cal_truth_metrics import cal_truth_pairs
from scripts.F1_related import cal_F1_related
from IO.write_results import write_results
from IO.write_stats import write_stats
from IO.block_start_end import block_start_end
from process.prefilter import prefilter
from process.extract_target import extract_target



PWD = os.path.dirname(os.path.realpath(__file__))

@click.command()
@click.option("-i", "--input", required = True, help = "Input vcf/bcf file for evaluation")
@click.option("-n", "--name", default = None, help = "User defined sample name of the input vcf file, [default: `extract from vcf`]")
@click.option("-c", "--compare", required = True, help = "Truth vcf/bcf file for comparsion")
@click.option("-r", "--ref", required = True, help = "Reference file fasta or fasta.fai")
@click.option("-o", "--output", required = True, help = "Output tsv file prefix, path can be added before the prefix, such as -o /test/output_name")
@click.option("-t", "--threads", default = 24, help = "Maximum numbers of parallel threads")
@click.option("--bed", default = None, help = r"Regions to only include, defined in bed file")
@click.option("--block", default = None, help = r"Output phasing block start and end to bed file, such as --block blocks.bed")
@click.option("--min-sv", default = 30, help = "Minimal length of Structral Variant")
@click.option("--chrom", default = ','.join(["chr" + str(i) for i in range(1, 23)]), help = "Chromosome to evaluate,use comma to join chromosome name e.g. --chrom chr1,chr2,chr3 [default:chr1-chr22,]")
@click.option("--sexchrom", default = ["chrX, chrY"], help = "Sex chromosme,use comma to join chromosome name e.g. --sexchrom chrX,chrY [default:chrX,chrY,]")
@click.option("--mincount", default = 2, help = "Minimum numbers of phased sites in a phase block [default: 2]")
@click.option("--no-sex", is_flag = True, help = "Ignore sex chromosome")
@click.option("--canonical", is_flag = True, help = "Canonical mode, only evaluate single mutation SNV")
@click.option("--only-snv", is_flag = True, help = "Only evaluate Single Nucleotide Variation")
@click.option("--no-sv", is_flag = True, help = "Ignore Structural Variant")
@click.option("--no-indel", is_flag = True, help = "Ignore insertion and deletion")
@click.option("--no-double", is_flag = True, help = "Ignore double heterozygous site")
@click.option("--verbose", is_flag = True, help = "Verbose mode print intermediate results to stdout")
@click.version_option(version="es-0.5.0", prog_name = r"phasing all-in-one evaluator(pie), based on Python 3.7+")
def main(input, name, compare, ref, output, threads, bed, block, min_sv, chrom, sexchrom, mincount, canonical, no_sex, only_snv, no_sv, no_indel, no_double, verbose):
    start_time = time.time()
    
    logging.basicConfig(level = logging.DEBUG, format = "%(asctime)s - %(levelname)s - %(message)s")
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    else:
        logging.getLogger().setLevel(logging.WARNING)

    # get aboslute path
    logging.info(f"Processing input parameters")
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
        clean_chrom = []
        for i in chrom:
            if i not in sexchrom:
                clean_chrom.append(i)
        chrom = clean_chrom[:]
    
    # sort chromosome
    chrom = sort_chrom(chrom)

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
    if not bed:
        target = None
    else:
        logging.info("Benchmarking of specified bed area")
        target = read_bed(bed, chrom)
        target_len = len(target)
        if target_len <= 0:
            logging.error("No region specified in the BED file. Please ensure the BED file follows the standard tab-separated format")
            sys.exit(1)
        else:
            logging.info(f"{target_len} region(s) specified in bed file")
    
    # read VCF files
    logging.info("Reading query VCF file")
    with Pool(threads) as p:
        query_vcf = p.starmap(read_vcf, [(input, c, target, min_sv) for c in chrom])

    logging.info("Reading truth VCF file")
    with Pool(threads) as p:
        truth_vcf = p.starmap(read_vcf, [(compare, c, target, min_sv) for c in chrom])


    logging.info("Evalating genotype and intersecting blocks")
    with Pool(threads) as p:
        intersect_result = p.starmap(intersect, [(query_vcf[i], truth_vcf[i], chrom[i]) for i in range(len(chrom))])
    blocks = [i[0] for i in intersect_result]
    genotype_result = [i[1] for i in intersect_result]

    loging.info("Writing genotype result")
    write_genotype(genotype_result, output, chrom)   

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
    
    
    # Output to files
    logging.info("Writing output file")
    if (not name):
        name = get_sample_name(input)
    write_results(evaluation_results, chrom, (NG50, NG90), output, name, verbose)
    

    logging.info("ALL DONE")   
    end_time = time.time()
    logging.info(f"Total processing time is {end_time - start_time} seconds.")




if __name__ == "__main__":
    main()

