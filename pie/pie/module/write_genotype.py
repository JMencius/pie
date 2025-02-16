import os
import sys
import math
import logging
from pie.module.F1_related import cal_all



def write_genotype(genotype: dict, prefix: str, chrom: list, bed_target: list):
    if not bed_target:
        with open(os.path.abspath(prefix + ".variant.stats.csv"), 'w') as f:
            header = ["Chromosome", "TP", "FP", "FN", "Precision", "Recall", "F1-score", "Heterozygous phased count", "Heterozygous phased percentage"]
            f.write(','.join(header))
            f.write('\n')
        
            TP_sum, FP_sum, FN_sum, phased_sum, unphased_sum = 0, 0, 0, 0, 0
            for i in range(len(chrom)):
                working_chrom = chrom[i]
                cdict = genotype[i]
                TP_sum += cdict["TP"]
                FP_sum += cdict["FP"]
                FN_sum += cdict["FN"]
                phased_sum += cdict["PC"]
                unphased_sum += cdict["UC"]

                cprecision, crecall, cf1 = cal_all(cdict["TP"], cdict["FP"], cdict["FN"])
            
                if (cdict["PC"] + cdict["UC"]) != 0:
                    phased_percentage = cdict["PC"] / (cdict["PC"] + cdict["UC"])
                else:
                    phased_percentage = math.nan

                to_write = [working_chrom, cdict["TP"], cdict["FP"], cdict["FN"], cprecision, crecall, cf1, cdict["PC"], phased_percentage]
                f.write(','.join([str(i) for i in to_write]))
                f.write('\n')
            total_precision, total_recall, total_f1 = cal_all(TP_sum, FP_sum, FN_sum)
            if (phased_sum + unphased_sum) != 0:
                total_phased_percentage = phased_sum / (phased_sum + unphased_sum)
            else:
                logging.error("No phased heterozygous variant found")
                total_phased_percentage = math.nan

            to_write = ["ALL", TP_sum, FP_sum, FN_sum, total_precision, total_recall, total_f1, phased_sum, total_phased_percentage]
            f.write(','.join([str(i) for i in to_write]))
            f.write('\n')
    else:
        with open(os.path.abspath(prefix + ".variant.stats.csv"), 'w') as f:
            header = ["Chromosome", "Start", "End", "TP", "FP", "FN", "Precision", "Recall", "F1-score", "Heterozygous phased count", "Heterozygous phased percentage"]
            f.write(','.join(header))
            f.write('\n')
            
            TP_sum, FP_sum, FN_sum, phased_sum, unphased_sum = 0, 0, 0, 0, 0
            for i in range(len(bed_target)):
                working_chrom = bed_target[i][0]
                cdict = genotype[i]
                TP_sum += cdict["TP"]
                FP_sum += cdict["FP"]
                FN_sum += cdict["FN"]
                phased_sum += cdict["PC"]
                unphased_sum += cdict["UC"]

                cprecision, crecall, cf1 = cal_all(cdict["TP"], cdict["FP"], cdict["FN"])

                if (cdict["PC"] + cdict["UC"]) != 0:
                    phased_percentage = cdict["PC"] / (cdict["PC"] + cdict["UC"])
                else:
                    phased_percentage = math.nan

                to_write = [working_chrom, bed_target[i][1], bed_target[i][2], cdict["TP"], cdict["FP"], cdict["FN"], cprecision, crecall, cf1, cdict["PC"], phased_percentage]
                f.write(','.join([str(i) for i in to_write]))
                f.write('\n')

            total_precision, total_recall, total_f1 = cal_all(TP_sum, FP_sum, FN_sum)
            if (phased_sum + unphased_sum) != 0:
                total_phased_percentage = phased_sum / (phased_sum + unphased_sum)
            else:
                logging.error("No phased heterozygous variant found")
                total_phased_percentage = math.nan

            to_write = ["ALL", '', '', TP_sum, FP_sum, FN_sum, total_precision, total_recall, total_f1, phased_sum, total_phased_percentage]
            f.write(','.join([str(i) for i in to_write]))
            f.write('\n')


