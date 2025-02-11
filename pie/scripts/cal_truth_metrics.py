
def cal_pairs(total: int) -> int:
    if total < 1:
        return 0
    else:
        total = int(total)
        return int((total * (total - 1)) / 2)


def cal_truth_pairs(truth_blocks: dict) -> int:
    pair_count = 0
    for block in truth_blocks:
        phased_count = 0
        for record in truth_blocks[block]:
            if record.isphased:
                phased_count += 1

        block_pairs = cal_pairs(phased_count)
        pair_count += block_pairs
    
    return pair_count

