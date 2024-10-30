
def cal_pairs(total: int) -> int:
    if total < 1:
        return 0
    else:
        total = int(total)
        return int((total * (total - 1)) / 2)


def cal_truth_pairs(truth_blocks: dict) -> int:
    pair_count = 0
    for block in truth_blocks:
        block_pairs = cal_pairs(len(truth_blocks[block]))
        pair_count += block_pairs
    
    return pair_count

