from scripts.myrecord import myrecord
from scripts.myblock import myblock
import Levenshtein
from collections import deque



def preprocess_truth(truth_dict: dict, mincount: int) -> dict:
    preprocess = dict()
    block_index = 0
    for phaseblock in truth_dict:
        if len(truth_dict[phaseblock]) >= mincount:
            count = 0
            for r in truth_dict[phaseblock]:
                if r.isphased:
                    preprocess[r.pos] = (r.ref, r.alt, str(block_index), r.left, r.right)
                    count += 1

        block_index += 1
    return preprocess



def clean_blocks(indict: dict, mincount: int) -> list:
    clean_list = list()
    for i in indict.values():
        if i.count >= mincount:
            i.length = i.end - i.start + 1
            clean_list.append(i)
    return clean_list


def find_start_end_idx(block_phasestat) -> tuple:
    if sum(block_phasestat) == 0:
        return None

    

def intersect(query: dict, truth: dict, chrom: str, mincount: int, min_sv: int) -> list:
    blocks = list()

    pre_truth = preprocess_truth(truth, mincount)
    

    for phaseblock in query:
        if len(query[phaseblock]) < mincount:
            continue

        block_dict = dict()
        for r in query[phaseblock]:
            if (r.pos in pre_truth):
                in_truth = pre_truth[r.pos]
                if r.ref == in_truth[0] and r.alt == in_truth[1]:
                    
                    block_index = in_truth[2]
                    # add block to hash table
                    if block_index not in block_dict:
                        block_dict[block_index] = myblock(chrom, block_index)
                    
                    # add things to block
                    block_dict[block_index].phasestat.append(int(r.isphased))
                    block_dict[block_index].variants_pos.append(r.pos)
                    block_dict[block_index].left.append(r.left)
                    block_dict[block_index].right.append(r.right)
                    block_dict[block_index].truthleft.append(in_truth[3])
                    block_dict[block_index].truthright.append(in_truth[4])
                    block_dict[block_index].weight.append(1)
                    
                    # add according to category
                    if r.isphased:
                        block_dict[block_index].phasecount += 1
                        if r.category == "SNV":
                            block_dict[block_index].snv += 1
                        else:
                            if r.category == "INDEL":
                                block_dict[block_index].indel += 1
                            else:
                                if r.category == "SV":
                                    block_dict[block_index].sv += 1

    
        cleaned = clean_blocks(block_dict, mincount)
        for i in cleaned:
            blocks.append(i)

    # sort according to start position
    blocks.sort(key = lambda K : K.start)

    return blocks



