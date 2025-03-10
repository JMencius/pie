import os


def write_block(blocks: list, output: str, mincount: int) -> None:
    with open(os.path.abspath(output + ".blocks.bed"), 'w') as f:
        # write each block
        blocks_start_end = dict()
        for i in blocks:
            for j in i.values():
                if (j.snv + j.indel + j.sv) >= mincount:
                    if j.chrom not in blocks_start_end:
                        blocks_start_end[j.chrom] = list()
                    blocks_start_end[j.chrom].append((j.start, j.end))
        
        for b in blocks_start_end.values():
            b.sort(key = lambda K : K[0])


        for c in blocks_start_end:
            for b in blocks_start_end[c]:
                f.write('\t'.join([c, str(b[0]), str(b[1])]))
                f.write('\n')

