def cal_precision(TP: int, FP: int) -> float:
    precision = TP / (TP + FP)
    return precision


def cal_recall(TP: int, FN: int) -> float:
    recall = TP / (TP + FN)
    return recall


def cal_F1(precision: float, recall: float) -> float:
    F1 = (2 * precision * recall) / (precision + recall)
    return F1


def cal_F1_related(total_pairs: int, eva_result: dict) -> dict:
    count = 0
    for c in eva_result:
        FP = c['pairwise_switch_error']
        TP = c['pairwise_event'] - FP
        FN = total_pairs[count] - c['pairwise_event']
        precision = cal_precision(TP, FP)
        recall = cal_recall(TP, FN)
        F1 = cal_F1(precision, recall)
        
        c['TP'] = TP
        c['FP'] = FP
        c['FN'] = FN
        c['precision'] =  precision
        c['recall'] = recall
        c['F1'] = F1

    return eva_result
