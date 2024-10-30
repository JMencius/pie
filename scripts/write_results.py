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
    with open(os.path.abspath(prefix + ".perchrom.tsv"), 'w') as f:
        f.write("Chromosome\t\
                Total Phased\t\
                Total SNV\t\
                Total INDEL\t\
                Total SV\t\
                Block count\t\
                Phased block length median\t\
                Total switch error\t\
                Switch error rate\t\
                Pairwise switch error\t\
                Pairwise event\t\
                Pairwise switch error rate\t\
                Pairwise precision\t\
                Pairwise recall\t\
                Pairwise F1\t\
                Total GHD\t\
                Total GHD event\t\
                GHD percentage\n")
        for c in range(len(eva_list)):
            f.write(f"chr{chrom[c]}\t\
                    {eva_list[c]['total_phase']}\t\
                    {eva_list[c]['total_snv']}\t\
                    {eva_list[c]['total_indel']}\t\
                    {eva_list[c]['total_sv']}\t\
                    {eva_list[c]['total_block']}\t\
                    {eva_list[c]['length_median']}\t\
                    {eva_list[c]['total_se']}\t\
                    {eva_list[c]['total_se'] / eva_list[c]['total_phase']}\t\
                    {eva_list[c]['pairwise_switch_error']}\t\
                    {eva_list[c]['pairwise_event']}\t\
                    {eva_list[c]['pairwise_switch_error'] / (eva_list[c]['pairwise_event'] - eva_list[c]['total_block'])}\t\
                    {eva_list[c]['precision']}\t\
                    {eva_list[c]['recall']}\t\
                    {eva_list[c]['F1']}\t\
                    {eva_list[c]['total_present']}\t\
                    {eva_list[c]['total_ghd']}\t\
                    {eva_list[c]['total_present'] / eva_list[c]['total_ghd']}\n")
    
    out_dict = cal_overall(eva_list)
    # write overall results
    with open(os.path.abspath(prefix + ".overall.tsv"), 'w') as g:
        g.write(f"Sample name\t\
                Total Phased\t\
                Total SNV\t\
                Total INDEL\t\
                Total SV\t\
                Total block\t\
                NG50\t\
                NG90\t\
                Switch Error count\t\
                Switch Error rate\t\
                Pairwise switch error\t\
                Pairwise event\t\
                Pairwise switch error rate\t\
                Pairwise precision\t\
                Pairwise recall\t\
                Pairwise F1\t\
                Total GHD\t\
                Total GHD event\t\
                GHD percentage\n")
        g.write(f"{sample_name}\t\
                {out_dict['total_phase']}\t\
                {out_dict['overall_snv']}\t\
                {out_dict['overall_indel']}\t\
                {out_dict['overall_sv']}\t\
                {out_dict['overall_block_count']}\t\
                {NG[0]}\t\
                {NG[1]}\t\
                {out_dict['se_count']}\t\
                {out_dict['se_rate']}\t\
                {out_dict['total_pse']}\t\
                {out_dict['total_pse_event']}\t\
                {out_dict['total_pse_rate']}\t\
                {out_dict['total_precision']}\t\
                {out_dict['total_recall']}\t\
                {out_dict['total_F1']}\t\
                {out_dict['total_ghd']}\t\
                {out_dict['total_ghd_event']}\t\
                {out_dict['ghd_percentage']}\n")
    

