import pyfastx
import os
import sys


def get_ref_len(ref: str, chrom: list) -> dict:
    suffix = os.path.splitext(ref)[1]
    if suffix == ".fa" or suffix == ".fasta":
        length_dict = fasta_mode(ref, chrom)
        return length_dict
    elif suffix == ".fai":
        length_dict = fai_mode(ref, chrom)
        return length_dict
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


def fasta_mode(ref: str, chrom: list) -> dict:
    chr_len = dict()
    exist = set()
    for name, seq in pyfastx.Fasta(ref, build_index = False):
        clean_name = clean(name)
        if clean_name in chrom:
            chr_len[clean_name] = len(seq)
            exist.add(clean_name)

    if len(chr_len) == len(chrom):
        return chr_len
    else:
        for i in chrom:
            if i not in exist:
                print(f"CRITICAL ERROR chr{i} not in reference file")
        print("Some chrosome missing in reference file, NG50 and NG90 will not be calculated")
        return -1


def fai_mode(ref_fai: str, chrom: list) -> dict:
    chr_len = dict()
    exist = set()
    with open(ref_fai, 'r') as f:
        for line in f:
            m = line.split()
            if len(m) >= 2:
                clean_name = clean(m[0])
                if clean_name in chrom:
                    chr_len[clean_name] = int(m[1])
                    exist.add(clean_name)

    if len(chr_len) == len(chrom):
        return chr_len
    else:
        for i in chrom:
            if i not in exist:
                print(f"CRITICAL ERROR chr{i} not in reference file")
        print("Some chrosome missing in reference file, NG50 and NG90 will not be calculated")
        return -1
                




