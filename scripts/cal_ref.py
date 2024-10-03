import pyfastx
import os
import sys


def cal_total_ref(ref: str, chrom: list) -> int:
    suffix = os.path.splitext(ref)[1]
    if suffix == ".fa" or suffix == ".fasta":
        total_length = fasta_mode(ref, chrom)
        return total_length
    elif suffix == ".fai":
        total_length = fai_mode(ref, chrom)
        return total_length
    else:
        raise ValueError("Invalid file suffix of reference file")


def clean(in_name: str) -> str:
    out_name = ""
    if len(in_name) >= 4:
        if in_name[: 3] == "chr" or in_name[: 3] == "Chr":
            out_name = in_name[3 :]

    if out_name:
        return out_name
    else:
        return in_name


def fasta_mode(ref: str, chrom: list) -> int:
    chr_len = list()
    exist = set()
    for name, seq in pyfastx.Fasta(ref, build_index = False):
        clean_name = clean(name)
        if clean_name in chrom:
            chr_len.append(len(seq))
            exist.add(clean_name)

    if len(chr_len) == len(chrom):
        return sum(chr_len)
    else:
        for i in chrom:
            if i not in exist:
                print(f"CRITICAL WARNING chr{i} not in reference file")
        print("Some chrosome missing in reference file, NG50 and NG90 will not be calculated")
        return -1


def fai_mode(ref_fai: str, chrom: list) -> int:
    chr_len = list()
    exist = set()
    with open(ref_fai, 'r') as f:
        for line in f:
            m = line.split()
            if len(m) >= 2:
                clean_name = clean(m[0])
                if clean_name in chrom:
                    chr_len.append(int(m[1]))
                    exist.add(clean_name)

    if len(chr_len) == len(chrom):
        return sum(chr_len)
    else:
        for i in chrom:
            if i not in exist:
                print(f"CRITICAL WARNING chr{i} not in reference file")
        print("Some chrosome missing in reference file, NG50 and NG90 will not be calculated")
        return -1
                




