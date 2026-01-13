from pie.module.pie_class import block
from pie.module.F1_related import cal_all
import ctypes
import os
import sys
import logging
import itertools
import numpy as np
from numba import njit

# load C so

script_dir = os.path.dirname(os.path.realpath(__file__))
lib = ctypes.CDLL(f"{script_dir}/lib/hamming.so")

# hamming distance C function
lib.hamming_distance.restype = ctypes.c_int
lib.hamming_distance.argtypes = [ctypes.c_char_p, ctypes.c_char_p]
def c_hamming_distance(str1, str2) -> int:
    if len(str1) != len(str2):
        raise ValueError("Length of the two input strings should be equal.")
    return lib.hamming_distance(str1.encode(), str2.encode())


# hamming comparision C function
lib.hamming_comparison.restype = ctypes.POINTER(ctypes.c_int)
lib.hamming_comparison.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_int]
# free result C function
lib.free_result.argtypes = [ctypes.POINTER(ctypes.c_int)]
lib.free_result.restype = None
def hamming_comparison(str1, str2):
    if len(str1) != len(str2):
        raise ValueError("Length of the two input strings should be equal.")
    
    length = len(str1)
    
    result_ptr = lib.hamming_comparison(str1.encode(), str2.encode(), length)
    result = [result_ptr[i] for i in range(length)]
    lib.free_result(result_ptr)
    
    return result



def cal_se_debug(hamming_list: list, variant_pos: list, working_chr: str) -> tuple:
    """
    calculate traditional switch error(SE)
    """
    event_count = len(hamming_list) - 1
    # 'S' means Start
    flag = 'S'
    se_count = 0
    c = 0
    for i in hamming_list:
        if i == 1:
            if flag == 0:
                se_count += 1
                print(f"{working_chr} {variant_pos[c]}-{variant_pos[c - 1]} SE")
                flag = 1

            if flag == 'S':
                flag = 1


        if i == 0:
            if flag == 1:
                se_count += 1
                print(f"{working_chr} {variant_pos[c]}-{variant_pos[c - 1]} SE")
                flag = 0
            if flag == 'S':
                flag = 0
        c += 1

    return (se_count, event_count)


def cal_se(hamming_list: list) -> tuple:
    """
    calculate traditional switch error(SE)
    """
    event_count = len(hamming_list) - 1
    # 'S' means Start
    flag = 'S'
    se_count = 0
    for i in hamming_list:
        if i == 1:
            if flag == 0:
                se_count += 1
                flag = 1

            if flag == 'S':
                flag = 1


        if i == 0:
            if flag == 1:
                se_count += 1
                flag = 0
            if flag == 'S':
                flag = 0
    
    return (se_count, event_count)


@njit(fastmath=True, nogil=True)
def calc_internal_tp_fp(sites, v0, v1, max_len):
    tp = 0.0
    fp = 0.0
    n = len(sites)
    

    for i in range(n):
        site_i = sites[i]
        val0_i = v0[i]
        val1_i = v1[i]
        
        for j in range(i + 1, n):
            site_j = sites[j]
            if site_i > site_j:
                dist = site_i - site_j
            else:
                dist = site_j - site_i
                
            if dist <= max_len:
                w = 1.0
            else:
                w = 1.0 / (dist - max_len)

            val0_j = v0[j]
            val1_j = v1[j]
            
            if (val0_i * val0_j) == (val1_i * val1_j):
                tp += w
            else:
                fp += w
     
           
    return tp, fp



@njit(fastmath=True, nogil=True)
def fast_calc_weights(arr_a, arr_b, truth_a, truth_b, max_len):
    total_weight = 0.0
    n_a = len(arr_a)
    n_b = len(arr_b)

    for i in range(n_a):
        val_a = arr_a[i]
        t_a = truth_a[i]
        
        for j in range(n_b):
            val_b = arr_b[j]
            t_b = truth_b[j]
            
            if t_a != t_b or val_a == val_b:
                continue

            if val_a > val_b:
                dist = val_a - val_b
            else:
                dist = val_b - val_a

            if dist <= max_len:
                total_weight += 1.0
            else:
                total_weight += 1.0 / (dist - max_len)

    return total_weight

def evaluate_interblock(listA, listB, truth_dict, max_len):
    arr_a = np.array(listA, dtype=np.int64)
    arr_b = np.array(listB, dtype=np.int64)

    unique_labels = sorted(list(set(truth_dict.values())))
    
    label_map = {label: idx for idx, label in enumerate(unique_labels)}

    truth_a_list = [label_map[truth_dict[x]] for x in listA]
    truth_b_list = [label_map[truth_dict[x]] for x in listB]

    truth_a = np.array(truth_a_list, dtype=np.int8)
    truth_b = np.array(truth_b_list, dtype=np.int8)
    
    return fast_calc_weights(arr_a, arr_b, truth_a, truth_b, max_len)


def blockwise_evaluate(chrom_block: dict, ref_len_dict: dict, mincount: int, truth_dict: dict, max_len: int, target_bed_list: list, genotype_FP: list, lmdb: bool) -> dict:
    len_list = list()
    total_snv, total_indel, total_sv, total_phase = 0, 0, 0, 0
    total_block = 0
    SE_denom, SE = 0, 0
    HD_denom, HD = 0, 0
    PSE_denom, PSE = 0, 0
    lmdb_list = list()

    more_than_2variant = list()
    idx = 0
    blocks = list(chrom_block.values())
    # calculate traditional metrics
    for b in blocks:
        # filter out very small block
        if (b.snv + b.indel + b.sv) < mincount:
            idx += 1
            continue
        
        more_than_2variant.append(idx)

        # b for block
        total_block += 1
        present_chrom = b.chrom
        total_snv += b.snv
        total_indel += b.indel
        total_sv += b.sv
        total_phase += b.snv + b.indel + b.sv
        len_list.append(b.end - b.start + 1)
        
        # calculate hamming distance
        hd = c_hamming_distance(b.queryleft, b.truthleft)
        
        if hd > 0.5 * len(b.queryleft):
            hd = c_hamming_distance(b.queryleft, b.truthright)
            compare_list = hamming_comparison(b.queryleft, b.truthright)
            compare_subject = b.truthright
        else:
            compare_list = hamming_comparison(b.queryleft, b.truthleft)
            compare_subject = b.truthleft
        
        lmdb_list.append(compare_list)
        
        HD_denom += len(b.queryleft)
        HD += hd
        
        # calculate switch error rate
        se_count, se_event_count = cal_se(compare_list)
        ## se_count, se_event_count = cal_se_debug(compare_list, list(b.subject.keys()), b.chrom)
        SE += se_count
        SE_denom += se_event_count
    
        idx += 1

    
    pairwise_TP, pairwise_FP, pairwise_FN = 0, 0, 0
    # calculate pairwise metrics
    newfn = set()
    for i in more_than_2variant:
        if len(blocks[i].querysymbol) >= 2 and len(blocks[i].truthsymbol) >= 2:
            v0_arr = np.array(blocks[i].querysymbol, dtype = np.int8)
            v1_arr = np.array(blocks[i].truthsymbol, dtype = np.int8)
            tp, fp = calc_internal_tp_fp(np.array(blocks[i].phase_variants, dtype = np.int64), v0_arr, v1_arr, max_len)
            pairwise_TP += tp
            pairwise_FP += fp

        fn1 = evaluate_interblock(list(blocks[i].unphase_variants), list(blocks[i].phase_variants), truth_dict, max_len)
        fn2 = evaluate_interblock(list(blocks[i].unphase_variants), list(blocks[i].unphase_variants), truth_dict, max_len)
        pairwise_FN += fn1 + fn2 / 2

    blocks_variants = list()
    merged = list()
    for i in blocks:
        temp = list(i.phase_variants) + list(i.unphase_variants)
        if len(temp) == 1:
            merged.extend(temp)
        else:
            blocks_variants.append(temp)

    blocks_variants.append(merged)

    for i in range(len(blocks_variants)):
        for j in range(i + 1, len(blocks_variants)):
            fn = evaluate_interblock(blocks_variants[i], blocks_variants[j], truth_dict, max_len)
            pairwise_FN += fn


    pairwise_precision, pairwise_recall, pairwise_f1 = cal_all(pairwise_TP, pairwise_FP, pairwise_FN)

    results={"total_phase": total_phase,
            "SE": SE,
            "SE_denom": SE_denom,
            "HD": HD,
            "HD_denom": HD_denom,
            "length_list": len_list,
            "snv": total_snv,
            "indel": total_indel,
            "sv": total_sv,
            "block_count": total_block,
            "PSE": PSE,
            "PSE_denom": PSE_denom,
            "pairwise_TP": pairwise_TP,
            "pairwise_FP": pairwise_FP,
            "pairwise_FN": pairwise_FN,
            "pairwise_precision": pairwise_precision,
            "pairwise_recall": pairwise_recall,
            "pairwise_f1": pairwise_f1,
            }
    

    if lmdb:
        lmdb_dict = dict()
        count = 0
        for k, b in chrom_block.items():
            # filter out very small block
            if (b.snv + b.indel + b.sv) < mincount:
                for pv in b.phase_variants:
                    lmdb_dict[pv] = {"ps": "unphase"}
                for uv in b.unphase_variants:
                    lmdb_dict[uv] = {"ps": "unphase"}
                continue    
            
            subject_hamming = lmdb_list[count]
            
   
            for uv in b.unphase_variants:
                lmdb_dict[uv] = {"ps": "unphase"}
            
            for pv, h in zip(b.phase_variants, subject_hamming):
                lmdb_dict[pv] = {"ps": k, "hc": h}
    
            count += 1


    if lmdb:
        return (results, lmdb_dict)
    else:
        return results



