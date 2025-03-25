import os


def write_block(blocks_start_end: list, output: str, chrom: list) -> None:
    with open(os.path.abspath(output + ".blocks.bed"), 'w') as f:
        # write each block
        for i in range(len(chrom)):
            c = chrom[i]
            for j in blocks_start_end[i]:
                f.write('\t'.join([c, str(j[0]), str(j[1])]))
                f.write('\n')

