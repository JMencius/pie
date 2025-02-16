import math


def cal_precision(TP: int, FP: int) -> float:
    denominator = TP + FP
    
    if denominator != 0:
        precision = TP / denominator
        return precision
    else:
        return math.nan


def cal_recall(TP: int, FN: int) -> float:
    denominator = TP + FN

    if denominator != 0:
        recall = TP / denominator
        return recall
    else:
        return math.nan


def cal_f1(precision: float, recall: float) -> float:
    denominator = precision + recall
    if (denominator != 0) and not(math.isnan(denominator)):
        F1 = (2 * precision * recall) / denominator
        return F1
    else:
        return math.nan


def cal_all(TP: int, FP: int, FN: int) -> tuple:
    precision = cal_precision(TP, FP)
    recall = cal_recall(TP, FN)
    f1 = cal_f1(precision, recall)

    return (precision, recall, f1)

