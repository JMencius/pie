class myblock:
    def __init__(self, chrom, idx):
        self.chrom: str = chrom
        self.idx: str = idx
        self.left: list = []
        self.right: list = []
        self.truthleft: list = []
        self.truthright: list = []
        self.weight: list = []
        self.length: int = 0
        self.count: int = 0
        self.snv: int = 0
        self.indel: int = 0
        self.sv: int = 0


    def __str__(self) -> str:
        return f"chromosome: {self.chrom} index: {self.idx} count: {self.count} length: {self.length} SNV: {self.snv} INDEL: {self.indel} SV: {self.sv}"


