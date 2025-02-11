import os

def write_stats(stat: list, prefix: str, chrom: list):
    with open(os.path.abspath(prefix + ".variant.stats.csv"), 'w') as f:
        header = ["Chromosome", "TP", "FP", "FN", "Precision", "Recall", "F1-score", "Heterozygous phased count", "Heterozygous phased percentage"]
        f.write(','.join(header))
        f.write('\n')
        
        TP_sum, FP_sum, FN_sum, phased_sum = 0, 0, 0, 0
        for i in range(len(chrom)):
            working_chrom = chrom[i]
            cdict = stat[i]
            TP_sum += cdict["TP"]
            FP_sum += cdict["FP"]
            FN_sum += cdict["FN"]
            phased_sum += cdict["TP_phased"]

            cprecision = cdict["TP"] / (cdict["TP"] + cdict["FP"])
            crecall = cdict["TP"] / (cdict["TP"] + cdict["FN"])
            cf1 = 2 * cprecision * crecall / (cprecision + crecall)
            to_write = [working_chrom, cdict["TP"], cdict["FP"], cdict["FN"], cprecision, crecall, cf1, cdict["TP_phased"], cdict["TP_phased"] / cdict["TP"]]
            f.write(','.join([str(i) for i in to_write]))
            f.write('\n')
        
        total_precision = TP_sum / (TP_sum + FP_sum)
        total_recall = TP_sum / (TP_sum + FN_sum)
        total_f1 = 2 * total_precision * total_recall / (total_precision + total_recall)
        to_write = ["ALL", TP_sum, FP_sum, FN_sum, total_precision, total_recall, total_f1, phased_sum, phased_sum / TP_sum]
        f.write(','.join([str(i) for i in to_write]))
        f.write('\n')


