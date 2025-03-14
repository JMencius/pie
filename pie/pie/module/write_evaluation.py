import os
import math
from pie.module.F1_related import cal_all
from pie.module.safediv import safediv



def write_evaluation(output: str, evaluation_results: list, chrom: list, name: str, NG50: int, NG90: int, bed: list) -> None:
    # write result in normal mode
    if not bed:
        with open(os.path.abspath(output + ".perchrom.csv"), 'w') as f:
            # write header
            perchrom_header = ["Chromosome", "Total phased", "SNV", "INDEL", "SV", "Block count", "NG50", "Total Switch error", "Switch error rate", "Hamming distance", "Hamming distance percentage", "Pairwise event", "Pairwise switch error rate", "Pairwise precision", "Pairwise recall", "Pairwise F1"]
            f.write(','.join(perchrom_header))
            f.write('\n')
            
            total_phase, total_snv, total_indel, total_sv, total_block = 0, 0, 0, 0, 0
            total_se, total_se_denom = 0, 0
            total_hd, total_hd_denom = 0, 0
            total_pse, total_pse_denom = 0, 0
            total_pairwise_TP, total_pairwise_FP, total_pairwise_FN = 0, 0, 0

            # write each chromosome and summarize
            for i in range(len(chrom)):
                s = evaluation_results[i]
                result = [chrom[i], s["total_phase"], s["snv"], s["indel"], s["sv"], s["block_count"], s["NG50"], s["SE"], safediv(s["SE"], s["SE_denom"]), s["HD"], safediv(s["HD"], s["HD_denom"]), s["PSE_denom"], safediv(s["PSE"], s["PSE_denom"]), s["pairwise_precision"], s["pairwise_recall"], s["pairwise_f1"]]
                f.write(','.join([str(t) for t in result]))
                f.write('\n')

                total_phase += s["total_phase"]
                total_snv += s["snv"]
                total_indel += s["indel"]
                total_sv += s["sv"]
                total_block += s["block_count"]

                total_se += s["SE"]
                total_se_denom += s["SE_denom"]

                total_hd += s["HD"]
                total_hd_denom += s["HD_denom"]

                total_pse += s["PSE"]
                total_pse_denom += s["PSE_denom"]

                total_pairwise_TP += s["pairwise_TP"]
                total_pairwise_FP += s["pairwise_FP"]
                total_pairwise_FN += s["pairwise_FN"]

        with open(os.path.abspath(output + ".overall.csv"), 'w') as f:
            # write header
            overall_header = ["Sample", "Total phased", "SNV", "INDEL", "SV", "Block count", "NG50", "NG90", "Total Switch error", "Switch error rate", "Hamming distance", "Hamming distance percentage", "Pairwise event", "Pairwise switch error rate", "Pairwise precision", "Pairwise recall", "Pairwise F1"]
            f.write(','.join(overall_header))
            f.write('\n')

            # calculate and write metrics
            overall_precision, overall_recall, overall_f1 = cal_all(total_pairwise_TP, total_pairwise_FP, total_pairwise_FN)
            overall_result = [name, total_phase, total_snv, total_indel, total_sv, total_block, NG50, NG90, total_se, safediv(total_se, total_se_denom), total_hd, safediv(total_hd, total_hd_denom), total_pse_denom, safediv(total_pse, total_pse_denom), overall_precision, overall_recall, overall_f1]
            f.write(','.join([str(t) for t in overall_result]))
            f.write('\n')
    
    # write result in bed mode
    else:
        with open(os.path.abspath(output + ".perregion.csv"), 'w') as f:
            # write header
            perchrom_header = ["Chromosome", "Start", "End", "Total phased", "SNV", "INDEL", "SV", "Block count", "Total Switch error", "Switch error rate", "Hamming distance", "Hamming distance percentage", "Pairwise event", "Pairwise switch error rate", "Pairwise precision", "Pairwise recall", "Pairwise F1"]
            f.write(','.join(perchrom_header))
            f.write('\n')
            
            total_phase, total_snv, total_indel, total_sv, total_block = 0, 0, 0, 0, 0
            total_se, total_se_denom = 0, 0
            total_hd, total_hd_denom = 0, 0
            total_pse, total_pse_denom = 0, 0
            total_pairwise_TP, total_pairwise_FP, total_pairwise_FN = 0, 0, 0

            # write each chromosome and summarize
            for i in range(len(bed)):
                s = evaluation_results[i]
                result = [bed[i][0], bed[i][1], bed[i][2], s["total_phase"], s["snv"], s["indel"], s["sv"], s["block_count"], s["SE"], safediv(s["SE"], s["SE_denom"]), s["HD"], safediv(s["HD"], s["HD_denom"]), s["PSE_denom"], safediv(s["PSE"], s["PSE_denom"]), s["pairwise_precision"], s["pairwise_recall"], s["pairwise_f1"]]
                f.write(','.join([str(t) for t in result]))
                f.write('\n')

                total_phase += s["total_phase"]
                total_snv += s["snv"]
                total_indel += s["indel"]
                total_sv += s["sv"]
                total_block += s["block_count"]

                total_se += s["SE"]
                total_se_denom += s["SE_denom"]

                total_hd += s["HD"]
                total_hd_denom += s["HD_denom"]

                total_pse += s["PSE"]
                total_pse_denom += s["PSE_denom"]

                total_pairwise_TP += s["pairwise_TP"]
                total_pairwise_FP += s["pairwise_FP"]
                total_pairwise_FN += s["pairwise_FN"]

        with open(os.path.abspath(output + ".overall.csv"), 'w') as f:
            # write header
            overall_header = ["Sample", "Total phased", "SNV", "INDEL", "SV", "Block count", "Total Switch error", "Switch error rate", "Hamming distance", "Hamming distance percentage", "Pairwise event", "Pairwise switch error rate", "Pairwise precision", "Pairwise recall", "Pairwise F1"]
            f.write(','.join(overall_header))
            f.write('\n')

            # calculate and write metrics
            overall_precision, overall_recall, overall_f1 = cal_all(total_pairwise_TP, total_pairwise_FP, total_pairwise_FN)
            overall_result = [name, total_phase, total_snv, total_indel, total_sv, total_block, total_se, safediv(total_se, total_se_denom), total_hd, safediv(total_hd, total_hd_denom), total_pse_denom, overall_precision, overall_recall, overall_f1]
            f.write(','.join([str(t) for t in overall_result]))
            f.write('\n')

















