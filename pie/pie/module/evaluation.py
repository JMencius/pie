from pie.module.pie_class import block
from pie.module.F1_related import cal_all
import ctypes
import os
import sys
import logging




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



def cal_pse(query: str, truth: str, phase_variants, max_len: int) -> tuple:
    """
    calculate pairwise swtich error(pse)
    """
    q, t = query, truth
    pse = 0
    event = 0
    cut_point = 1
    current_pos = 0
    while current_pos < len(q) - 1:
        while (cut_point < len(q) - 1) and ((phase_variants[cut_point] - phase_variants[current_pos]) <= max_len):
            cut_point += 1

        q1 = q[current_pos]
        t1 = t[current_pos]
        
        subject_q = q[current_pos + 1 : cut_point + 1]
        subject_t = t[current_pos + 1 : cut_point + 1]

        dist = c_hamming_distance(subject_q, subject_t)
        #print(dist)
        if q1 != t1:
            dist = len(subject_q) - dist
        event += len(subject_q)
        pse += dist
        
        current_pos += 1

    return (pse, event)



def Cn2(n: int) -> int:
    return n * (n - 1) / 2



def process_truth_count(truth_count: dict, max_len: int) -> int:
    total_pairs = 0
    for pos_list in truth_count.values():
        if len(pos_list) >= 2:
            pos_list.sort()
            j = 1
            for i in range(len(pos_list) - 1):
                while (j < len(pos_list) - 1) and ((pos_list[j] - pos_list[i]) <= max_len):
                    j += 1
                total_pairs += (j - i)


    return total_pairs


def blockwise_evaluate(chrom_block: dict, ref_len_dict: dict, mincount: int, truth_count: dict, max_len: int, target_bed_list: list, present_chrom: str) -> dict:
    len_list = list()
    total_snv, total_indel, total_sv, total_phase = 0, 0, 0, 0
    total_block = 0
    SE_denom, SE = 0, 0
    HD_denom, HD = 0, 0
    PSE_denom, PSE = 0, 0
    pairwise_FN = 0

    # calculate total pairs
    total_pairs = process_truth_count(truth_count, max_len)

    for b in chrom_block.values():
        # filter out very small block
        if (b.snv + b.indel + b.sv) < mincount:
            continue

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

        HD_denom += len(b.queryleft)
        HD += hd
        
        # calculate switch error rate
        se_count, se_event_count = cal_se(compare_list)
        SE += se_count
        SE_denom += se_event_count
        
        # calculate pairwise switch error
        pse_count, pse_event_count = cal_pse(b.queryleft, compare_subject, b.phase_variants, max_len)
        PSE += pse_count
        PSE_denom += pse_event_count
        

    pairwise_FP = PSE
    pairwise_TP = PSE_denom - PSE
    pairwise_FN = total_pairs - PSE_denom


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
    

    return results



