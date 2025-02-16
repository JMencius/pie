import os


def write_block(blocks: list, output: str, chrom: list) -> None:
    with open(os.path.abspath(output + ".blocks.bed"), 'w') as f:
        # write header
        f.write('\t'.join(["Chromosome", "Start", "End"]))
        f.write('\n')

        # write each block
        for i in range(len(chrom)):
            c = chrom[i]
            blocks = []
            for j in blocks[i].values():
                blocks.append((j.start, j.end))

            blocks.sort(key = lambda K : K[0])

            for b in blocks:
                f.write('\t'.join([c, str(b[0]), str(b[1])]))
                f.write('\n')

