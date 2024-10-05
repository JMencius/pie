from scripts.myrecord import myrecord
from scripts.myblock import myblock
import Levenshtein



def preprocess_truth(truth_dict: dict) -> dict:
    preprocess = dict()
    block_index = 0
    for phaseblock in truth_dict:
        for r in truth_dict[phaseblock]:
            preprocess[r.pos] = (r.ref, r.alt, str(block_index), r.left, r.right)

        block_index += 1
    return preprocess



def clean_blocks(inlist: list, mincount: int) -> list:
    clean_list = list()
    for i in inlist:
        if i.count > mincount:
            clean_list.append(i)
    return clean_list



def intersect(query: dict, truth: dict, chrom: str, mincount: int, min_sv: int) -> list:
    blocks = list()

    pre_truth = preprocess_truth(truth)

    for phaseblock in query:
        if len(query[phaseblock]) < mincount:
            continue

        tempblock = myblock(chrom, "")
        minmax = [None, None]
        for r in query[phaseblock]:
            if r.pos in pre_truth:
                in_truth = pre_truth[r.pos]
                if r.ref == in_truth[0] and r.alt == in_truth[1]:
                    if tempblock.idx == "":
                        tempblock.idx = in_truth[2]
                    else:
                        if in_truth[2] != tempblock.idx:
                            blocks.append(tempblock)
                            tempblock = myblock(chrom, in_truth[2])
                    tempblock.left.append(r.left)
                    tempblock.right.append(r.right)
                    tempblock.truthleft.append(in_truth[3])
                    tempblock.truthright.append(in_truth[4])
                    tempblock.count += 1
                    if (not minmax[0]) and (not minmax[1]):
                        minmax[0] = r.pos
                        minmax[1] = r.pos
                    else:
                        minmax[0] = min(r.pos, minmax[0])
                        minmax[1] = max(r.pos, minmax[1])

                    if r.category == "SNV":
                        tempblock.weight.append(1)
                        tempblock.snv += 1
                    else:
                        alt_list = list(r.alt)
                        if len(r.alt) == 1:
                            weight = Levenshtein.distance(r.ref, str(alt_list[0]))
                        else:
                            weigth = Levenshtein.distance(str(alt_list[0]), str(alt_list[1]))
                        
                        tempblock.weight.append(weight)
                        if r.category == "INDEL":
                            tempblock.indel += 1
                        else:
                            if r.category == "SV":
                                tempblock.sv += 1
        if minmax[0] and minmax[1]:
            tempblock.length = minmax[1] - minmax[0] + 1
            blocks.append(tempblock)

    cleaned = clean_blocks(blocks, mincount)
    return cleaned



