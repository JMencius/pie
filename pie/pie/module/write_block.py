import os


def write_block(blocks_start_end: list, output: str, chrom: list, bed_target_list) -> None:
    if not bed_target_list:
        with open(os.path.abspath(output + ".blocks.bed"), 'w') as f:
            # write each block
            for i in range(len(chrom)):
                c = chrom[i]
                for j in blocks_start_end[i]:
                    f.write('\t'.join([c, str(j[0]), str(j[1])]))
                    f.write('\n')
    else:
        with open(os.path.abspath(output + ".blocks.bed"), 'w') as f:
            # write each block
            for i in blocks_start_end:
                f.write('\t'.join([str(i[0]), str(i[1]), str(i[2])]))
                f.write('\n')

