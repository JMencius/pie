def cal_NGx0(len_list: list, total_len: dict, p: int) -> int:
    total_len = sum(total_len.values())
    target = total_len * p / 100

    if sum(len_list) < target:
        return None
    
    cummulative = 0
    for i in len_list:
        cummulative += i
        if cummulative >= target:
            return i
