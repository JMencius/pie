import os
from scripts.F1_related import cal_precision
from scripts.F1_related import cal_recall
from scripts.F1_related import cal_F1


def cal_overall(eva_list: list) -> dict:
    out_dict = dict()
    total_snv, total_indel, total_sv = 0, 0, 0
    total_block = 0
    total_event, total_se = 0, 0
    total_ghd_event, total_ghd = 0, 0
    total_pse, total_pse_event = 0, 0
    total_TP, total_FP, total_FN = 0, 0, 0
    for i in eva_list:
        total_snv += i['total_snv']
        total_indel += i['total_indel']
        total_sv += i['total_sv']
        total_block += i['total_block']
        total_event += i['total_phase']
        total_se += i['total_se']
        total_ghd_event += i['total_ghd']
        total_ghd += i['total_present']
        total_pse += i['pairwise_switch_error']
        total_pse_event += i['pairwise_event']
        total_TP += i['TP']
        total_FP += i['FP']
        total_FN += i['FN']
    
    total_precision = cal_precision(total_TP, total_FP)
    total_recall = cal_recall(total_TP, total_FN)
    total_F1 = cal_F1(total_precision, total_recall)

    out_dict["total_phase"] = total_event
    out_dict["overall_snv"] = total_snv
    out_dict["overall_indel"] = total_indel
    out_dict["overall_sv"] = total_sv
    out_dict["overall_block_count"] = total_block
    out_dict["se_count"] = total_se
    out_dict["se_rate"] = total_se / (total_event - total_block)
    out_dict["total_ghd"] = total_ghd
    out_dict["total_ghd_event"] = total_ghd_event
    out_dict["ghd_percentage"] = total_ghd / total_ghd_event
    out_dict["total_pse"] = total_pse
    out_dict["total_pse_event"] = total_pse_event
    out_dict["total_pse_rate"] = total_pse / total_pse_event
    out_dict["total_TP"] = total_TP
    out_dict["total_FP"] = total_FP
    out_dict["total_FN"] = total_FN
    out_dict["total_precision"] = total_precision
    out_dict["total_recall"] = total_recall
    out_dict["total_F1"] = total_F1

    return out_dict


def write_results(eva_list: list, chrom: list, NG: tuple, prefix: str, sample_name: str, verbose: bool):
    # write per chromosome results
    with open(os.path.abspath(prefix + ".perchrom.csv"), 'w') as f:
        header = ["Chromosome", "Total Phased", "Total SNV", "Total INDEL", "Total SV", "Block count", "Phased block length median", "Total switch error", "Switch error rate", "Pairwise switch error", "Pairwise event", "Pairwise switch error rate", "Pairwise precision", "Pairwise recall", "Pairwise F1", "Hamming distance", "Hamming distance percentage"]
        f.write(','.join(header))
        f.write('\n')

        for c in range(len(eva_list)):
            temp = [f"chr{chrom[c]}", eva_list[c]['total_phase'], eva_list[c]['total_snv'], eva_list[c]['total_indel'], eva_list[c]['total_sv'], eva_list[c]['total_block'], eva_list[c]['NG50'], eva_list[c]['total_se'], eva_list[c]['total_se'] / eva_list[c]['total_phase'], eva_list[c]['pairwise_switch_error'], eva_list[c]['pairwise_event'], eva_list[c]['pairwise_switch_error'] / (eva_list[c]['pairwise_event'] - eva_list[c]['total_block']), eva_list[c]['precision'], eva_list[c]['recall'], eva_list[c]['F1'], eva_list[c]['total_present'], eva_list[c]['total_present'] / eva_list[c]['total_ghd']]
            f.write(','.join([str(i) for i in temp]))
            f.write('\n')
    
    out_dict = cal_overall(eva_list)
    # write overall results
    with open(os.path.abspath(prefix + ".overall.csv"), 'w') as g:
        header = ["Sample name", "Total Phased", "Total SNV", "Total INDEL", "Total SV", "Total block", "NG50", "NG90", "Switch Error count", "Switch Error rate", "Pairwise switch error", "Pairwise event", "Pairwise switch error rate", "Pairwise precision", "Pairwise recall", "Pairwise F1", "Hamming distance", "Hamming distance percentage"]
        g.write(','.join(header))
        g.write('\n')

        result = [sample_name, out_dict['total_phase'], out_dict['overall_snv'], out_dict['overall_indel'], out_dict['overall_sv'], out_dict['overall_block_count'], NG[0], NG[1], out_dict['se_count'], out_dict['se_rate'], out_dict['total_pse'], out_dict['total_pse_event'], out_dict['total_pse_rate'], out_dict['total_precision'], out_dict['total_recall'], out_dict['total_F1'], out_dict['total_ghd'], out_dict['ghd_percentage']]
        g.write(','.join([str(i) for i in result]))
        g.write('\n')
