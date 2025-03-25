from multiprocessing import Pool
from pie.module.cal_NGx0 import cal_NGx0
import logging



def get_raw_block_length(variant_dict: dict) -> list:
    raw_blocks = dict()
    for i in variant_dict.values():
        if i.isphased:
            if i.ps not in raw_blocks:
                raw_blocks[i.ps] = list()
            raw_blocks[i.ps].append(i.pos)
    
    block_len = list()
    block_start_end = list()
    for i in raw_blocks.values():
        block_start_end.append((min(i), max(i)))
        block_len.append(max(i) - min(i) + 1)
    
    # sort according to start position
    block_start_end.sort(key = lambda K: K[0])

    return (block_start_end, block_len)       



def cal_NG50(variants: list, len_dict: dict, chrom: list, threads: int) -> tuple:
    with Pool(threads) as p:
        temp = p.map(get_raw_block_length, [i for i in variants])
    
    block_start_end = [i[0] for i in temp]
    block_len = [i[1] for i in temp]
    

    NG50 = dict()
    total_chrom_len = 0
    total_block_len = list()
    for i in range(len(chrom)):
        working_chr = chrom[i]
        chrom_len = len_dict[working_chr]
        current_NG50 = cal_NGx0(block_len[i], chrom_len, 50)
        NG50[working_chr] = current_NG50
        
        total_chrom_len += chrom_len
        total_block_len = total_block_len + block_len[i]

    NG50["ALL"] = cal_NGx0(total_block_len, total_chrom_len, 50)


    return (block_start_end, NG50)
