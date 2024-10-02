import pyfastx
import os


def cal_half_ref(ref: str) -> float:
    suffix = os.path.splitext(ref)[1]
    if suffix == ".fa" or suffix == ".fasta":
        pass
    elif suffix == ".fai":
        pass
    else:
        raise ValueError("Invalid file suffix of reference file")
