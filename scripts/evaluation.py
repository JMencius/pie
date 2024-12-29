from scripts.overall_metrics import cal_NGx0
import ctypes
import os
import sys
from scripts.myblock import myblock
from collections import deque
from typing import Deque


##N50, NG50, NG90, Switch Error, Generalize hamming distance, phasing percentage



# load C so

script_dir = os.path.dirname(os.path.realpath(__file__))
lib = ctypes.CDLL(f"{script_dir}/hamming.so")

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



def cal_SE(hamming_list: list) -> tuple:
    """
    calculate traditional switch error(SE)
    """
    event_count = len(hamming_list)
    flag = 'S'
    se_count = 0
    if event_count == 0:
        raise ValueError("Length of the evaluating block is 0")
    else:
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


def cal_pse(query: Deque[str], truth: Deque[str]) -> tuple:
    """
    calculate pairwise swtich error(pse)
    """
    q = query.copy()
    t = truth.copy()
    error = 0
    event = 0
    while len(q) > 1:
        #print(q, t)
        q_left = q.popleft()
        t_left = t.popleft()
        dist = c_hamming_distance(''.join(q), ''.join(t))
        #print(dist)
        if q_left != t_left:
            dist = len(q) - dist
        event += len(q)
        error += dist

    return (error, event)


def cal_GHD(weight_list: list, hamming_list: list) -> tuple:
    """
    calculate generalize hamming distance(GHD)
    """
    total_ghd = sum(weight_list)
    present_ghd = 0
    for i in range(len(weight_list)):
        present_ghd += weight_list[i] * hamming_list[i]
    
    return (present_ghd, total_ghd)



def blockwise_evaluate(chrom_block: list, idx: int, ref_len_dict: dict, verbose: bool) -> dict:
    # return tuple is (median, Switch Error, Generalized Hamming Distance)
    result_list = list()
    len_list = list()
    total_snv, total_indel, total_sv = 0, 0, 0
    total_block = 0
    total_phase = 0
    total_se, total_event = 0, 0
    total_present, total_ghd = 0, 0
    pse, pse_event = 0, 0
    present_chrom = None
    for in_block in chrom_block:
        present_chrom = in_block.chrom
        total_block += 1
        total_snv += in_block.snv
        total_indel += in_block.indel
        total_sv += in_block.sv
        

        total_phase += in_block.count
        len_list.append(in_block.length)
        query: str = ''.join(in_block.left)
        truth: str = ''.join(in_block.truthleft)
    
        hd = c_hamming_distance(query, truth)

        if hd > 0.5 * len(query):
            query: str = ''.join(in_block.right)
            hd = c_hamming_distance(query, truth)
        compare_list = hamming_comparison(query, truth)

        # calculate switch error rate
        se_count, event_count = cal_SE(compare_list)
        total_se += se_count
        total_event += event_count

        # calculate pairwise switch error
        a, b = cal_pse(in_block.left, in_block.truthleft)
        pse += a
        pse_event += b
        
        # calculate hamming distance
        current_ghd, current_present = cal_GHD(in_block.weight, compare_list)
        total_present += current_ghd
        total_ghd += current_present
        

    if len(len_list) > 0:
        NG50 = cal_NGx0(len_list, ref_len_dict[present_chrom], 50)
    else:
        length_median = None
        print(f"CRITIAL WARNING: No phase block in chromosome {present_chrom}, please check the vcfs")

    return {"total_phase": total_phase,
            "NG50": NG50,
            "total_se": total_se,
            "total_event": total_event,
            "total_present": total_present,
            "total_ghd": total_ghd,
            "length_list": len_list,
            "total_snv": total_snv,
            "total_indel": total_indel,
            "total_sv": total_sv,
            "total_block": total_block,
            "pairwise_switch_error": pse,
            "pairwise_event": pse_event}




"""
if __name__ == "__main__":
    str1 = "1@001"
    str2 = "10011"
    distance = c_hamming_distance(str1, str2)
    print(f"Hamming distance: {distance}")
"""


