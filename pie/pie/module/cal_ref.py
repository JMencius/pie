import pyfastx
import os
import sys
import logging


def get_ref_len(ref: str, chrom: list) -> dict:
    suffix = os.path.splitext(ref)[1]
    if suffix == ".fa" or suffix == ".fasta":
        length_dict = fasta_mode(ref, chrom)
        return length_dict
    elif suffix == ".fai":
        length_dict = fai_mode(ref, chrom)
        return length_dict
    else:
        logging.error("Invalid file suffix of the reference file, must be either .fasta/.fa or .fai")
        raise ValueError("Invalid file suffix of the reference file")



def fasta_mode(ref: str, chrom: list) -> dict:
    chr_len = dict()
    exist = set()
    for name, seq in pyfastx.Fasta(ref, build_index = False):
        if name in chrom:
            chr_len[name] = len(seq)
            exist.add(name)

    if len(chr_len) == len(chrom):
        return chr_len
    else:
        for i in chrom:
            if i not in exist:
                logging.warning(f"{i} not in reference file")
        return chr_len



def fai_mode(ref_fai: str, chrom: list) -> dict:
    chr_len = dict()
    exist = set()
    with open(ref_fai, 'r') as f:
        for line in f:
            m = line.split()
            if len(m) >= 2:
                name = m[0]
                if name in chrom:
                    chr_len[name] = int(m[1])
                    exist.add(name)

    if len(chr_len) == len(chrom):
        return chr_len
    else:
        for i in chrom:
            if i not in exist:
                logging.warning(f"{i} not in reference file")
        return chr_len
                




