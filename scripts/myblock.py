class myblock:
    def __init__(self, chrom, idx, left, right, truthleft, truthright, weight):
        self.chrom: str = chrom
        self.idx: str = idx
        self.left: list = left
        self.right: list = right
        self.truthleft: list = truthleft
        self.truthright: list = truthright
        self.weight: list = weight
        self.count: int = 0
