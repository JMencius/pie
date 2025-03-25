import logging
import sys


def read_bed(bedfile: str, chrom: list) -> tuple:
    out_dict = dict()
    region_count = 0
    with open(bedfile, 'r') as f:
        for line in f:
            m = (line.strip()).split('\t')
            if len(m) >= 3:
                if m[0] in chrom:
                    if m[0] not in out_dict:
                        out_dict[m[0]] = list()
                    if m[1].isdigit() and m[2].isdigit():
                        region_count += 1
                        out_dict[m[0]].append((int(m[1]), int(m[2])))

    sorted_out_dict = dict()
    for i, j in out_dict.items():
        temp = sorted(j, key = lambda K : K[0])
        sorted_out_dict[i] = temp
    
    if target_len <= 0:
        logging.error("No region specified in the BED file. Please ensure the BED file follows the standard tab-separated format")
        sys.exit(1)
    else:
        logging.info(f"{target_len} region(s) specified in bed file")

    return (region_count, sorted_out_dict)
