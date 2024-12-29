def read_bed(bedfile: str) -> list:
    out_list = []
    with open(bedfile, 'r') as f:
        for line in f:
            m = (line.strip()).split('\t')
            if len(m) >= 3:
                if ("chr" in m[0]) or ("Chr" in m[0]):
                    out_list.append([m[0][3: ], int(m[1]), int(m[2])])
                else:
                    out_list.append([m[0], int(m[1]), int(m[2])])

    return out_list
