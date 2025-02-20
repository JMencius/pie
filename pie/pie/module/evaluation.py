from pie.module.pie_class import block
from pie.module.overall_metrics import cal_NGx0
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



def cal_pse(query: str, truth: str) -> tuple:
    """
    calculate pairwise swtich error(pse)
    """
    q, t = query, truth
    pse = 0
    event = 0
    while len(q) > 1:
        #print(q, t)
        q1 = q[0]
        t1 = t[0]
        
        q = q[1:]
        t = t[1:]

        dist = c_hamming_distance(q, t)
        #print(dist)
        if q1 != t1:
            dist = len(q) - dist
        event += len(q)
        pse += dist

    return (pse, event)



def cal_pairwise_FN(blocks: dict, mincount: int) -> int:
    break_cause = 0

    start = True
    last_block_variants = 0
    total_variants = 0
    total_genotypeFN = 0
    for b in blocks.values():
        vc = b.snv + b.indel + b.sv
        # filter out very small block
        if vc < mincount:
            continue

        total_variants += vc
        total_genotypeFN += b.FN
        if start:
            last_block_variants = vc
            start = False
        else:
            break_cause += (vc * total_variants) * FN_correct_coefficient(vc, total_variants)

    genotypeFN_cause = total_variants * total_genotypeFN + (total_genotypeFN - 1) * total_genotypeFN / 2

    return genotypeFN_cause + break_cause


def FN_correct_coefficient(n1: int, n2: int) -> float:
    correct_coefficient = 4 / (n1 + n2 - 2) * (1 + (min(n1, n2) - (n1 + n2) // 2) / sum([n1, n2]))
    return correct_coefficient



def blockwise_evaluate(chrom_block: dict, ref_len_dict: dict, mincount: int, target_bed_list: list) -> dict:
    len_list = list()
    total_snv, total_indel, total_sv, total_phase = 0, 0, 0, 0
    total_block = 0
    SE_denom, SE = 0, 0
    HD_denom, HD = 0, 0
    PSE_denom, PSE = 0, 0
    present_chrom = None
    pairwise_FN = 0
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
        pse_count, pse_event_count = cal_pse(b.queryleft, compare_subject)
        PSE += pse_count
        PSE_denom += pse_event_count
        
        # calculate pairwise recall
        unphase_count = b.FN
        if unphase_count != 0:
            pairwise_FN += (b.snv + b.indel + b.sv) * unphase_count + ((unphase_count - 1) * unphase_count / 2)
    

    pairwise_FP = PSE
    pairwise_TP = PSE_denom - PSE
    pairwise_FN = cal_pairwise_FN(chrom_block, mincount)

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
            "pairwise_f1": pairwise_f1}
    
    if not target_bed_list:
        if len(len_list) > 0:
            NG50 = cal_NGx0(len_list, ref_len_dict[present_chrom], 50)
            results["NG50"] = NG50
        else:
            logging.warning(f"No phase block in chromosome {present_chrom}")
    

    return results



if __name__ == "__main__":
    str1 = "1@001"
    str2 = "10011"
    distance = c_hamming_distance(str1, str2)
    print(f"Hamming distance: {distance}")

    distance_list = hamming_comparison(str1, str2)
    print(cal_se(distance_list))


