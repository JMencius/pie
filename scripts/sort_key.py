def sort_key(s: str) -> int:
    flag = 1
    for i in s:
        if i.isalpha():
            flag = 0

    if flag == 1:
        return int(s)
    else:
        return ord(s[0])
