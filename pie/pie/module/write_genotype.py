import os
import sys
import math
import logging
from pie.module.F1_related import cal_all
from pie.module.safediv import safediv


def write_genotype(genotype: dict, phase_data: dict, prefix: str, chrom: list, bed_target: list):
    if not bed_target:
        with open(os.path.abspath(prefix + ".variant.stats.csv"), 'w') as f:
            header = ["Chromosome", "TP", "FP", "FN", "Precision", "Recall", "F1-score", "SNV phased count", "SNV total count", "SNV phased percentage", "INDEL phased count", "INDEL count", "INDEL phased percentage", "SV phased count", "SV total count", "SV phased percentage", "Total phased count", "Total count", "Total phased percentage",]
            f.write(','.join(header))
            f.write('\n')
        
            TP_sum, FP_sum, FN_sum, phased_sum, unphased_sum = 0, 0, 0, 0, 0
            snv_phased_count, snv_total_count = 0, 0
            indel_phased_count, indel_total_count = 0, 0
            sv_phased_count, sv_total_count = 0, 0
            for i in range(len(chrom)):
                working_chrom = chrom[i]
                cdict = genotype[i]
                pd_dict = phase_data[i]
                TP_sum += cdict["TP"]
                FP_sum += cdict["FP"]
                FN_sum += cdict["FN"]
                phased_sum += cdict["PC"]
                unphased_sum += cdict["UC"]

                snv_phased_count += pd_dict["SNV"][0]
                snv_total_count += (pd_dict["SNV"][0] + pd_dict["SNV"][1])

                indel_phased_count += pd_dict["INDEL"][0]
                indel_total_count += (pd_dict["INDEL"][0] + pd_dict["INDEL"][1])

                sv_phased_count += pd_dict["SV"][0]
                sv_total_count += (pd_dict["SV"][0] + pd_dict["SV"][1])

                cprecision, crecall, cf1 = cal_all(cdict["TP"], cdict["FP"], cdict["FN"])
            
                phased_percentage = safediv(cdict["PC"], (cdict["PC"] + cdict["UC"]))
                
                snv_phased_percentage = safediv(pd_dict["SNV"][0], (pd_dict["SNV"][0] + pd_dict["SNV"][1]))
                indel_phased_percentage = safediv(pd_dict["INDEL"][0], (pd_dict["INDEL"][0] + pd_dict["INDEL"][1]))
                sv_phased_percentage = safediv(pd_dict["SV"][0], (pd_dict["SV"][0] + pd_dict["SV"][1]))                

                to_write = [working_chrom, cdict["TP"], cdict["FP"], cdict["FN"], cprecision, crecall, cf1, \
                            pd_dict["SNV"][0], pd_dict["SNV"][0] + pd_dict["SNV"][1], snv_phased_percentage, \
                            pd_dict["INDEL"][0], pd_dict["INDEL"][0] + pd_dict["INDEL"][1], indel_phased_percentage, \
                            pd_dict["SV"][0], pd_dict["SV"][0] + pd_dict["SV"][1], sv_phased_percentage, \
                            cdict["PC"], cdict["PC"] + cdict["UC"],  phased_percentage]
                f.write(','.join([str(i) for i in to_write]))
                f.write('\n')

            total_precision, total_recall, total_f1 = cal_all(TP_sum, FP_sum, FN_sum)
            to_write = ["ALL", TP_sum, FP_sum, FN_sum, total_precision, total_recall, total_f1, \
                        snv_phased_count, snv_total_count, safediv(snv_phased_count, snv_total_count), \
                        indel_phased_count, indel_total_count, safediv(indel_phased_count, indel_total_count), \
                        sv_phased_count, sv_total_count, safediv(sv_phased_count, sv_total_count), \
                        phased_sum, phased_sum + unphased_sum, safediv(phased_sum, phased_sum + unphased_sum)]
            f.write(','.join([str(i) for i in to_write]))
            f.write('\n')
    else:
        with open(os.path.abspath(prefix + ".variant.stats.csv"), 'w') as f:
            header = ["Chromosome", "Start", "End", "TP", "FP", "FN", "Precision", "Recall", "F1-score", "SNV phased count", "SNV total count", "SNV phased percentage", "INDEL phased count", "INDEL count", "INDEL phased percentage", "SV phased count", "SV total count", "SV phased percentage", "Total phased count", "Total count", "Total phased percentage",]
            f.write(','.join(header))
            f.write('\n')
            
            TP_sum, FP_sum, FN_sum, phased_sum, unphased_sum = 0, 0, 0, 0, 0
            snv_phased_count, snv_total_count = 0, 0
            indel_phased_count, indel_total_count = 0, 0
            sv_phased_count, sv_total_count = 0, 0

            for i in range(len(bed_target)):
                working_chrom = bed_target[i][0]
                cdict = genotype[i]
                pd_dict = phase_data[i]
                TP_sum += cdict["TP"]
                FP_sum += cdict["FP"]
                FN_sum += cdict["FN"]
                phased_sum += cdict["PC"]
                unphased_sum += cdict["UC"]
                
                snv_phased_count += pd_dict["SNV"][0]
                snv_total_count += (pd_dict["SNV"][0] + pd_dict["SNV"][1])

                indel_phased_count += pd_dict["INDEL"][0]
                indel_total_count += (pd_dict["INDEL"][0] + pd_dict["INDEL"][1])

                sv_phased_count += pd_dict["SV"][0]
                sv_total_count += (pd_dict["SV"][0] + pd_dict["SV"][1])

                cprecision, crecall, cf1 = cal_all(cdict["TP"], cdict["FP"], cdict["FN"])

                phased_percentage = safediv(cdict["PC"], (cdict["PC"] + cdict["UC"]))
                
                snv_phased_percentage = safediv(pd_dict["SNV"][0], (pd_dict["SNV"][0] + pd_dict["SNV"][1]))
                indel_phased_percentage = safediv(pd_dict["INDEL"][0], (pd_dict["INDEL"][0] + pd_dict["INDEL"][1]))
                sv_phased_percentage = safediv(pd_dict["SV"][0], (pd_dict["SV"][0] + pd_dict["SV"][1]))
                
                to_write = [working_chrom, bed_target[i][1], bed_target[i][2], cdict["TP"], cdict["FP"], cdict["FN"], cprecision, crecall, cf1, \
                            pd_dict["SNV"][0], pd_dict["SNV"][0] + pd_dict["SNV"][1], snv_phased_percentage, \
                            pd_dict["INDEL"][0], pd_dict["INDEL"][0] + pd_dict["INDEL"][1], indel_phased_percentage, \
                            pd_dict["SV"][0], pd_dict["SV"][0] + pd_dict["SV"][1], sv_phased_percentage, \
                            cdict["PC"], cdict["PC"] + cdict["UC"],  phased_percentage]

                f.write(','.join([str(i) for i in to_write]))
                f.write('\n')

            total_precision, total_recall, total_f1 = cal_all(TP_sum, FP_sum, FN_sum)

            total_phased_percentage = safediv(phased_sum, (phased_sum + unphased_sum))

            to_write = ["ALL", '', '', TP_sum, FP_sum, FN_sum, total_precision, total_recall, total_f1, \
                        snv_phased_count, snv_total_count, safediv(snv_phased_count, snv_total_count), \
                        indel_phased_count, indel_total_count, safediv(indel_phased_count, indel_total_count), \
                        sv_phased_count, sv_total_count, safediv(sv_phased_count, sv_total_count), \
                        phased_sum, phased_sum + unphased_sum, safediv(phased_sum, phased_sum + unphased_sum)]
            f.write(','.join([str(i) for i in to_write]))
            f.write('\n')


