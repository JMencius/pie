import re

def remove_chr_prefix(s):
    return re.sub(r'^[cC]hr', '', s)

def custom_sort(s):
    if s.isdigit():
        return (0, int(s))
    else:
        return (1, s)

def sort_chrom(chrom: list) -> list:
    for_sort_list = [(remove_chr_prefix(i), i) for i in chrom]
    for_sort_list.sort(key = lambda K: custom_sort(K[0]))
    sorted_chrom = [i[1] for i in for_sort_list]

    return sorted_chrom

if __name__ == "__main__":
    print(sort_chrom(["chr1", "chr2", "chrX", "ChrY"]))
