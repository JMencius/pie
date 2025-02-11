def cal_NGx0(len_list: list, total_len: int, p: int) -> int:
    total_len = int(total_len)
    target = total_len * p / 100

    if sum(len_list) < target:
        return None
    
    cummulative = 0
    for i in len_list:
        cummulative += i
        if cummulative >= target:
            return i
