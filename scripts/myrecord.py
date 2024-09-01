class myrecord:
    def __init__(self, chrom, pos, ref, alt, ps, left, right, category):
        self.chrom : str = chrom
        self.pos : int = pos
        self.ref : str = ref
        self.alt : set = alt
        self.len : int = len(self.alt)
        self.ps : str = ps
        self.left : str = left
        self.right : str = right
        self.category : str = category

    def __str__(self) -> str:
        return (f"CHR:{self.chrom} POS:{self.pos} REF:{self.ref} ALT:{self.alt} LENGTH:{self.len} PS:{self.ps} LEFT: {self.left} RIGHT : {self.right} TYPE:{self.category}")

