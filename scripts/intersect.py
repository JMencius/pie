from scripts.myrecord import myrecord
from scripts.myblock import myblock
import Levenshtein



def preprocess_truth(truth_dict: dict, mincount: int) -> dict:
    preprocess = dict()
    block_index = 0
    for phaseblock in truth_dict:
        if len(truth_dict[phaseblock]) >= mincount:
            count = 0
            for r in truth_dict[phaseblock]:
                preprocess[r.pos] = (r.ref, r.alt, str(block_index), r.left, r.right)
                count += 1

        block_index += 1
    return preprocess



def clean_blocks(inlist: dict, mincount: int) -> list:
    clean_list = list()
    for i in inlist.values():
        if i.count >= mincount:
            i.length = i.end - i.start + 1
            clean_list.append(i)
    return clean_list



def intersect(query: dict, truth: dict, chrom: str, mincount: int, min_sv: int) -> list:
    blocks = list()

    pre_truth = preprocess_truth(truth, mincount)

    for phaseblock in query:
        if len(query[phaseblock]) < mincount:
            continue

        block_dict = dict()
        for r in query[phaseblock]:
            if r.pos in pre_truth:
                in_truth = pre_truth[r.pos]
                if r.ref == in_truth[0] and r.alt == in_truth[1]:
                    block_index = in_truth[2]
                    # add block to hash table
                    if block_index not in block_dict:
                        block_dict[block_index] = myblock(chrom, block_index)
                    
                    # add things to block
                    block_dict[block_index].left.append(r.left)
                    block_dict[block_index].right.append(r.right)
                    block_dict[block_index].truthleft.append(in_truth[3])
                    block_dict[block_index].truthright.append(in_truth[4])
                    block_dict[block_index].count += 1

                    # determine start end
                    if (not block_dict[block_index].start) and (not block_dict[block_index].end):
                        block_dict[block_index].start = r.pos
                        block_dict[block_index].end = r.pos
                    else:
                        block_dict[block_index].start = min(r.pos, block_dict[block_index].start)
                        block_dict[block_index].end = max(r.pos, block_dict[block_index].end)
                    
                    # add according to category
                    if r.category == "SNV":
                        block_dict[block_index].snv += 1
                        block_dict[block_index].weight.append(1)
                    else:
                        if r.category == "INDEL":
                            block_dict[block_index].indel += 1
                        else:
                            if r.category == "SV":
                                block_dict[block_index].sv += 1

                        alt_list = list(r.alt)
                        if len(r.alt) == 1:
                            weight = Levenshtein.distance(r.ref, str(alt_list[0]))
                        else:
                            weight = abs(Levenshtein.distance(str(alt_list[0]), str(alt_list[1])))
                        
                        block_dict[block_index].weight.append(weight)
    
        cleaned = clean_blocks(block_dict, mincount)
        for i in cleaned:
            blocks.append(i)

    # sort according to start position
    blocks.sort(key = lambda K : K.start)

    return blocks



