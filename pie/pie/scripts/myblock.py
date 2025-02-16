from collections import deque
from typing import Deque


class myblock:
    def __init__(self, chrom, idx):
        self.chrom: str = chrom
        self.idx: str = idx
        self.phasestat: Deque[int] = deque()
        self.varints_pos: Deque[bool] = deque()
        self.left:  Deque[str] = deque()
        self.right: Deque[str] = deque()
        self.truthleft: Deque[str] = deque()
        self.truthright: Deque[str] = deque()
        self.weight: Deque[int] = deque()
        self.length: int = 0
        self.phasecount: int = 0
        self.snv: int = 0
        self.indel: int = 0
        self.sv: int = 0
        self.start: int = None
        self.end: int = None


    def __str__(self) -> str:
        return f"chromosome: {self.chrom} Start: {self.start} End: {self.end} index: {self.idx} phasecount: {self.phasecount} length: {self.length} SNV: {self.snv} INDEL: {self.indel} SV: {self.sv}"


