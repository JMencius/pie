def read_bed(bedfile: str, chrom: list) -> dict:
    out_dict = dict()
    with open(bedfile, 'r') as f:
        for line in f:
            m = (line.strip()).split('\t')
            if len(m) >= 3:
                if m[0] in chrom:
                    if m[0] not in out_dict:
                        out_dict[m[0]] = list()
                    if m[1].isdigit() and m[2].isdigit():
                        out_dict[m[0]].append((int(m[1]), int(m[2])))

    sorted_out_dict = dict()
    for i, j in out_dict.items():
        temp = sorted(j, key = lambda K : K[0])
        sorted_out_dict[i] = temp

    return sorted_out_dict
