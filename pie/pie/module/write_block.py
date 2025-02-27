import os


def write_block(blocks: list, output: str, chrom: list, mincount: int) -> None:
    with open(os.path.abspath(output + ".blocks.bed"), 'w') as f:
        # write each block
        for i in range(len(chrom)):
            c = chrom[i]
            blocks_start_end = []
            for j in blocks[i].values():
                if (j.snv + j.indel + j.sv) >= mincount:
                    blocks_start_end.append((j.start, j.end))

            blocks_start_end.sort(key = lambda K : K[0])

            for b in blocks_start_end:
                f.write('\t'.join([c, str(b[0]), str(b[1])]))
                f.write('\n')

