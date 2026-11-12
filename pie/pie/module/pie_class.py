from sortedcontainers import SortedDict
from collections import deque
from typing import Deque
import hashlib
import json


def make_variant_key(pos, ref, alt):
    """Return a stable per-chromosome key, independent of ALT ordering."""
    catalogue = json.dumps([ref, sorted(alt)], separators=(",", ":"))
    return (pos, hashlib.sha256(catalogue.encode("utf-8")).hexdigest())



class pievariant:
    def __init__(self, chrom, pos, ref, alt, ps, left, right, type, isphased, isdouble, key):
        self.chrom : str = chrom
        self.pos : int = pos
        self.ref : str = ref
        self.alt : list = alt
        self.ps : str = ps
        self.left : str = left
        self.right : str = right
        self.type : str = type
        self.isphased : bool = isphased
        self.isdouble: bool = isdouble
        self.key : tuple = key

    def __str__(self) -> str:
        return (f"CHR:{self.chrom} POS:{self.pos} REF:{self.ref} ALT:{self.alt} PS:{self.ps} LEFT: {self.left} RIGHT: {self.right} TYPE:{self.type} PHASED?: {self.isphased} DOUBLE_MUT?: {self.isdouble}")


class block:
    def __init__(self, chrom, index):
        # init static properties
        self.chrom: str = chrom
        self.index: tuple = index
        # dynamic properties
        self.subject = SortedDict()
        self.phase_variants: Deque[tuple] = deque()
        self.unphase_variants: Deque[tuple] = deque()
        self.FN: int = 0
        self.snv: int = 0
        self.indel: int = 0
        self.sv: int = 0
        # properties to finalize
        self.start: int = None
        self.end: int = None
        self.queryleft: str = ""
        self.truthleft: str = ""
        self.truthright: str = ""
        self.querysymbol: list = list()
        self.truthsymbol: list = list()

    def add_phased_variant(self, query_variant, truth_variant):
        if truth_variant.pos == query_variant.pos:
            variant_key = query_variant.key
            self.phase_variants.append(variant_key)
            # Encode the actual alleles consistently across both VCF records.
            q = [query_variant.ref] + query_variant.alt
            t = [truth_variant.ref] + truth_variant.alt
            qleft, qright = q[query_variant.left], q[query_variant.right]
            tleft, tright = t[truth_variant.left], t[truth_variant.right]
            if qleft == qright or {qleft, qright} != {tleft, tright}:
                raise ValueError("Phased comparison requires matching heterozygous alleles")
            first_allele = min(qleft, qright)
            qbit = int(qleft != first_allele)
            tbit = int(tleft != first_allele)
            self.subject[variant_key] = (qbit, tbit, 1 - tbit)

            if tbit == 0:
                self.truthsymbol.append(-1)
            else:
                self.truthsymbol.append(1)

            if qbit == 0:
                self.querysymbol.append(-1)
            else:
                self.querysymbol.append(1)

            if query_variant.type == "SNV":
                self.snv += 1
            elif query_variant.type == "INDEL":
                self.indel += 1
            else:
                self.sv += 1
    
    def add_unphased_variant(self, variant_key):
        self.unphase_variants.append(variant_key)
        self.FN += 1

    def finalize(self):
        self.phase_variants = deque(self.subject.keys())
        self.start = min(key[0] for key in self.phase_variants)
        self.end = max(key[0] for key in self.phase_variants)
        self.querysymbol = [2 * item[0] - 1 for item in self.subject.values()]
        self.truthsymbol = [2 * item[1] - 1 for item in self.subject.values()]
        self.queryleft = ''.join([str(i[0]) for i in self.subject.values()])
        self.truthleft = ''.join([str(i[1]) for i in self.subject.values()])
        self.truthright = ''.join([str(i[2]) for i in self.subject.values()])
           
        



        
