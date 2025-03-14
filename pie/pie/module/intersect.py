from pie.module.pie_class import pievariant
from pie.module.pie_class import block


def check_genotype(query_variant, truth_variant) -> int:
    """
    return 0: not the same genotype
    return 1: exact the same genotype
    """
    q = [query_variant.ref] + query_variant.alt
    t = [truth_variant.ref] + truth_variant.alt
    q_left, q_right = q[query_variant.left], q[query_variant.right]
    t_left, t_right = t[truth_variant.left], t[truth_variant.right]
    

    if {q_left, q_right} == {t_left, t_right}:
        return 1
    else:
        return 0


def check_filters(filters: dict, pievariant) -> bool:
    if filters["only_snv"]:
        if pievariant.type != "SNV":
            return False
    
    if filters["only_indel"]:
        if pievariant.type != "INDEL":
            return False

    if filters["only_sv"]:
        if pievariant.type != "SV":
            return False

    if filters["no_snv"]:
        if pievariant.type == "SNV":
            return False
    
    if filters["no_indel"]:
        if pievariant.type == "INDEL":
            return False

    if filters["no_sv"]:
        if pievariant.type == "SV":
            return False

    if filters["no_double"]:
        if pievariant.isdouble:
            return False

    return True


def intersect(query: dict, truth: dict, chrom: str, filters: dict) -> dict:
    phaseblock = dict()
    genotype_TP, genotype_FP, genotype_FN = 0, 0, 0
    total_phase_count, total_unphase_count, interblock_unphase = 0, 0, 0
    unphase_variant = set()

    # [phase_count, unphase_count]
    phase_data = {"SNV": [0, 0], "INDEL": [0, 0], "SV": [0, 0]}

    # intersect blocks and calculate truth pairs
    truth_blocks = dict()

    for i in truth:
        if i not in query:
            genotype_FN += 1
        else:
            query_variant = query[i]
            truth_variant = truth[i]
            check_result = check_genotype(query_variant, truth_variant)

            if check_result == 0:
                genotype_FP += 1
            else:
                if truth_variant.isphased and truth_variant.ref == query_variant.ref:
                    if check_filters(filters, query_variant) and check_filters(filters, truth_variant):
                        genotype_TP += 1
                        if truth_variant.ps not in truth_blocks:
                            truth_blocks[truth_variant.ps] = list()
                        truth_blocks[truth_variant.ps].append(truth_variant.pos)
                        
                        if not query_variant.isphased:
                            total_unphase_count += 1
                            unphase_variant.add(i)
                            phase_data[query_variant.type][1] += 1
                        else:
                            total_phase_count += 1
                            phase_data[query_variant.type][0] += 1
                            if (truth_variant.ps, query_variant.ps) not in phaseblock:
                                phaseblock[(truth_variant.ps, query_variant.ps)] = block(query_variant.chrom, (truth_variant.ps, query_variant.ps))
                            phaseblock[(truth_variant.ps, query_variant.ps)].add_phased_variant(query_variant, truth_variant)
    

    # finalize blocks
    for b in phaseblock.values():
        b.finalize()
    
    # fill unphase variant back to block
    blocks_start_end = list()
    for b in phaseblock.values():
        blocks_start_end.append((b.start, b.end))
    

    for site in unphase_variant:
        in_block = find_closest_block(site, blocks_start_end)
        if not in_block:
            interblock_unphase += 1
        else:
            for b in phaseblock.values():
                if (b.start, b.end) == in_block:
                    b.add_unphased_variant(site)
    genotype_result = {"TP": genotype_TP, "FP": genotype_FP, "FN": genotype_FN, "PC": total_phase_count, "UC": total_unphase_count, "IUC": interblock_unphase}

    return (phaseblock, genotype_result, truth_blocks, phase_data)



def find_closest_block(site: int, intervals: list):
    temp = dict()
    for i in intervals:
        if i[0] < site < i[1]:
            temp[i] = min(site - i[0], i[1] - site)
    
    if not temp:
        return None
    else:
        temp_list = list(temp.items())
        temp_list.sort(key = lambda K : K[1])
        return temp_list[0][0]


    
            






