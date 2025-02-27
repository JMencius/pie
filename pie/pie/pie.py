import click
import sys
import os
import time
import logging
from multiprocessing import Pool
from pie.module.sort_chrom import sort_chrom
from pie.module.cal_ref import get_ref_len
from pie.module.read_bed import read_bed
from pie.module.read_vcf import read_vcf
from pie.module.intersect import intersect
from pie.module.write_genotype import write_genotype
from pie.module.check_filters import check_filters
from pie.module.evaluation import blockwise_evaluate
from pie.module.overall_metrics import cal_NGx0
from pie.module.write_block import write_block
from pie.module.write_evaluation import write_evaluation


PWD = os.path.dirname(os.path.realpath(__file__))

@click.command()
@click.option("-i", "--input", required = True, help = "Input vcf/bcf file for evaluation")
@click.option("-n", "--name", default = "Sample", help = "User defined sample name of the input vcf file, [default: Sample]")
@click.option("-c", "--compare", required = True, help = "Truth vcf/bcf file for comparsion")
@click.option("-r", "--ref", required = True, help = "Reference file fasta or fasta.fai")
@click.option("-o", "--output", required = True, help = "Output tsv file prefix, path can be added before the prefix, such as -o /test/output_name")
@click.option("-t", "--threads", default = 24, help = "Maximum numbers of parallel threads")
@click.option("-m", "--max-len", default = 250 * 10**3, help = "Maximum length between variants for pairwise calculation")
@click.option("--bed", default = None, help = r"Regions to only include, defined in bed file")
@click.option("--block", is_flag = True, help = r"Output phasing block start and end to bed file in output")
@click.option("--min-sv", default = 30, help = "Minimal length of Structral Variant")
@click.option("--chrom", default = ','.join(["chr" + str(i) for i in range(1, 23)]), help = "Chromosome to evaluate,use comma to join chromosome name e.g. --chrom chr1,chr2,chr3 [default:chr1-chr22,]")
@click.option("--sexchrom", default = ["chrX, chrY"], help = "Sex chromosme,use comma to join chromosome name e.g. --sexchrom chrX,chrY [default:chrX,chrY,]")
@click.option("--mincount", default = 2, help = "Minimum numbers of phased sites in a phase block [default: 2]")
@click.option("--no-sex", is_flag = True, help = "Ignore sex chromosome")
@click.option("--canonical", is_flag = True, help = "Canonical mode, only evaluate single mutation SNV")
@click.option("--only-snv", is_flag = True, help = "Only evaluate single nucleotide variation")
@click.option("--only-indel", is_flag = True, help = "Only evaluate insertion and deletion")
@click.option("--only-sv", is_flag = True, help = "Only evaluate structrual variant")
@click.option("--no-snv", is_flag = True, help = "Ignore single nucleotide variation")
@click.option("--no-indel", is_flag = True, help = "Ignore insertion and deletion")
@click.option("--no-sv", is_flag = True, help = "Ignore structural variant")
@click.option("--no-double", is_flag = True, help = "Ignore double heterozygous site")
@click.option("--verbose", is_flag = True, help = "Verbose mode print intermediate results to stdout")
@click.version_option(version="es-0.6.0-e", prog_name = r"phasing all-in-one evaluator(pie), based on Python 3.7+")
def main(input, name, compare, ref, output, threads, max_len, bed, block, min_sv, chrom, sexchrom, mincount, canonical, no_sex, only_snv, only_indel, only_sv, no_snv, no_indel, no_sv, no_double, verbose):
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
    chrom = [i for i in chrom.split(',')]
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
        bed_target = None
    else:
        logging.info("Benchmarking of specified bed area")
        target_len, bed_target = read_bed(bed, chrom)
        if target_len <= 0:
            logging.error("No region specified in the BED file. Please ensure the BED file follows the standard tab-separated format")
            sys.exit(1)
        else:
            logging.info(f"{target_len} region(s) specified in bed file")
    
    # read VCF files
    logging.info("Reading query VCF file")
    with Pool(threads) as p:
        query_vcf = p.starmap(read_vcf, [(input, c, bed_target, min_sv) for c in chrom])

    logging.info("Reading truth VCF file")
    with Pool(threads) as p:
        truth_vcf = p.starmap(read_vcf, [(compare, c, bed_target, min_sv) for c in chrom])


    logging.info("Evaluating genotype and intersecting blocks")
    # operate in normal mode
    filters = {"only_snv": only_snv, "only_indel": only_indel, "only_sv": only_sv, "no_snv": no_snv, "no_indel": no_indel, "no_sv": no_sv, "no_double": no_double}
    filter_result = check_filters(filters)
    if not filter_result:
        sys.exit(1)
    
    if not bed_target:
        target_bed_list = None
        with Pool(threads) as p:
            intersect_result = p.starmap(intersect, [(query_vcf[i], truth_vcf[i], chrom[i], filters) for i in range(len(chrom))])
        
    else:
        # flatten dict of list into dict
        query_pool = dict()
        truth_pool = dict()
        for i in query_vcf:
            for j in i:
                query_pool[j] = i[j]
        for i in truth_vcf:
            for j in i:
                truth_pool[j] = i[j]
        # operate in bed mode
        target_bed_list = list(query_pool.keys())
        with Pool(threads) as p:
            intersect_result = p.starmap(intersect, [(query_pool[i], truth_pool[i], i[0], filters) for i in target_bed_list])

    blocks = [i[0] for i in intersect_result]
    genotype_result = [i[1] for i in intersect_result]
    truth_count = [i[2] for i in intersect_result]
    

    logging.info("Writing genotype result")
    write_genotype(genotype_result, output, chrom, target_bed_list) 
    
    long_sum, short_sum = 0, 0
    logging.info("Evluating blocks")
    with Pool(threads) as r:
        evaluation_results = r.starmap(blockwise_evaluate, [(blocks[i], len_dict, mincount, truth_count[i], max_len, target_bed_list) for i in range(len(blocks))])

    if not target_bed_list:
        logging.info("Calculating overall NG50 and NG90")
        phase_len = list()
        for eva in evaluation_results:
            phase_len += eva["length_list"]
        phase_len.sort(reverse = True)
        NG50 = cal_NGx0(phase_len, total_ref_len, 50)
        NG90 = cal_NGx0(phase_len, total_ref_len, 90)
    else:
        NG50, NG90 = None, None

    ## print(evaluation_results[0], NG50, NG90)

    # Output to files
    logging.info("Writing output file")
    if block:
        ### need to filter
        logging.info("Writing block bed file")
        write_block(blocks, output, chrom, mincount)
    
    logging.info("Writing evaluation output file")
    write_evaluation(output, evaluation_results, chrom, name, NG50, NG90, target_bed_list)
    

    logging.info("ALL DONE")   
    end_time = time.time()
    logging.info(f"Total processing time is {end_time - start_time} seconds.")


if __name__ == "__main__":
    main()

