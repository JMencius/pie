import click
import sys
import os
import time
import logging
from multiprocessing import Pool
from pie.module.cli import cli
from pie.module.cal_ref import get_ref_len
from pie.module.read_bed import read_bed
from pie.module.read_vcf import read_vcf
from pie.module.intersect import intersect
from pie.module.write_genotype import write_genotype
from pie.module.check_filters import check_filters
from pie.module.evaluation import blockwise_evaluate
from pie.module.raw_metrics import cal_NG50
from pie.module.raw_metrics import get_raw_block_length
from pie.module.write_block import write_block
from pie.module.write_evaluation import write_evaluation
from pie.module.write_lmdb import write_lmdb



def main():
    start_time = time.time()
    
    param = cli(standalone_mode = False)     
   
    if isinstance(param, int):
        sys.exit(0)
    input, name, compare, ref, output, threads, max_len, bed, min_sv, chrom, sexchrom, mincount, canonical, block, no_sex, only_snv, only_indel, only_sv, no_snv, no_indel, no_sv, no_double, no_sort, include_genotype, lmdb, verbose = param

    logging.basicConfig(level = logging.DEBUG, format = "%(asctime)s - %(levelname)s - %(message)s")
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    else:
        logging.getLogger().setLevel(logging.WARNING)

    # read reference .fa or .fai
    logging.info("Reading reference file")
    len_dict = get_ref_len(ref, chrom)


    # get target area
    if not bed:
        bed_target = None
    else:
        logging.info("Reading specified bed areas")
        target_len, bed_target = read_bed(bed, chrom, no_sort)

    
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
    if filter_result == False:
        sys.exit(1)    
  
    if not bed_target:
        target_bed_list = None
        with Pool(threads) as p:
            intersect_result = p.starmap(intersect, [(query_vcf[i], truth_vcf[i], chrom[i], filters, include_genotype) for i in range(len(chrom))])
        
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
            intersect_result = p.starmap(intersect, [(query_pool[i], truth_pool[i], i[0], filters, include_genotype) for i in target_bed_list])

    blocks = [i[0] for i in intersect_result]
    genotype_result = [i[1] for i in intersect_result]
    truth_count = [i[2] for i in intersect_result]
    phase_count = [i[3] for i in intersect_result]   
    genotype_FP = [i[4] for i in intersect_result]


    logging.info("Writing genotype result")
    write_genotype(genotype_result, phase_count, output, chrom, target_bed_list) 

    
    logging.info("Evluating phase blocks")
    with Pool(threads) as r:
        evaluation_results = r.starmap(blockwise_evaluate, [(blocks[i], len_dict, mincount, truth_count[i], max_len, target_bed_list, include_genotype, genotype_FP[i], lmdb) for i in range(len(blocks))])

    if lmdb:
        lmdb_list = [i[1] for i in evaluation_results]
        evaluation_results = [i[0] for i in evaluation_results]

        logging.info("Outputting intermediate result to lmdb")
        write_lmdb(lmdb_list, chrom, output)              


    logging.info("Calculating NG50 metrics")
    if not bed_target:
        raw_blocks_start_end, NG50 = cal_NG50(query_vcf, len_dict, chrom, threads)
    else:
        raw_blocks_start_end = list()
        for i,j in query_pool.items():
            tem, _ = get_raw_block_length(j)        
            
            for t in tem:
                raw_blocks_start_end.append([i[0], t[0], t[1]])
        
        NG50 = None
    # Output to files
    logging.info("Writing output file")
    if block:
        logging.info("Writing block bed file")
        write_block(raw_blocks_start_end, output, chrom, target_bed_list)
    

    logging.info("Writing evaluation output file")
    write_evaluation(output, evaluation_results, chrom, name, NG50, target_bed_list)
    

    logging.info("ALL DONE")   
    end_time = time.time()
    logging.info(f"Total processing time is {end_time - start_time} seconds.")


if __name__ == "__main__":
    main()


