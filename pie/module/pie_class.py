from sortedcontainers import SortedDict
from collections import deque
from typing import Deque


class variant:
    def __init__(self, chrom, pos, ref, alt, ps, left, right, type, isphased, isdouble):
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

    def __str__(self) -> str:
        return (f"CHR:{self.chrom} POS:{self.pos} REF:{self.ref} ALT:{self.alt} PS:{self.ps} LEFT: {self.left} RIGHT: {self.right} TYPE:{self.type} PHASED?: {self.isphased} DOUBLE_MUT?: {self.isdouble}")


class block:
    def __init__(self, chrom, index):
        # init static properties
        self.chrom: str = chrom
        self.index: tuple = index
        # dynamic properties
        self.subject = SortedDict()
        self.phase_variants: Deque[int] = deque()
        self.unphase_variants: Deque[int] = deque()
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

    def add_phased_variant(self, query_variant, truth_variant, action: int):
        if truth_variant.pos == query_variant.pos:
            self.phase_variants.append(query_variant.pos)
            if action == 1:
                self.subject[query_variant.pos] = (query_variant.left, truth_variant.left, truth_variant.right)
            elif action == -1:
                self.subject[query_variant.pos] = (query_variant.right, truth_variant.left, truth_variant.right)

            if query_variant.type == "SNV":
                self.snv += 1
            elif query_variant.type == "INDEL":
                self.indel += 1
            else:
                self.sv += 1
    
    def add_unphased_variant(self, unphase_pos):
        self.unphase_variants.append(unphase_pos)
        self.FN += 1

    def finalize(self):
        self.start = phase_variants[0]
        self.end = phase_variants[-1]
        self.queryleft = ''.join([i[0] for i in self.subject.values()])
        self.truthleft = ''.join([i[1] for i in self.subject.values()])
        self.truthright = ''.join([i[2] for i in self.subject.values()])

           
        



        
