from scripts.myrecord import myrecord


def convert(instr: str):
    if ',' in instr:
        return set(instr[1].split(','))
    else:
        return instr


def get_category(inlist: list, min_sv: int) -> tuple:
    double = 0
    if ',' in inlist[1]:
        double = 1
        alt_list = inlist[1].split(',')
    else:
        alt_list = [inlist[1]]

    alt_list = [inlist[0]] + alt_list
    left = alt_list[int(inlist[2])]
    right = alt_list[int(inlist[3])]
    max_len = max(len(left), len(right))

    if len(inlist[0]) == 1:
        if max_len == 1:
            return (double, "SNV")
        elif 1 <= max_len <= min_sv:
            return (double, "INDEL")
        else:
            return (double, "SV")
    elif 1 < len(inlist[0]) <= min_sv:
        if 1 <= max_len <= min_sv:
            return (double, "INDEL")
        else:
            return (double, "SV")
    else:
        return (double, "SV")


def get_flag(isdouble: int, category: str, only_snv: bool, no_sv: bool, no_indel: bool, no_double: bool) -> bool:
    if no_double and isdouble == 1:
        return False

    if only_snv and category != "SNV":
        return False

    if no_sv and category == "SV":
        return False

    if no_indel and category == "INDEL":
        return False

    return True


#                            0    1    2        3         4        5
#  variants[int(var_pos)] = (ref, alt, gt_left, gt_right, isphase, ps)

def prefilter(query: dict, truth: dict, working_chr: str, min_sv: int, only_snv: bool, no_sv: bool, no_indel: bool, no_double: bool):
    TP, FP, FN = 0, 0, 0
    TP_phased = 0
    query_blocks = dict()
    truth_blocks = dict()

    for site in truth:
        if site not in query:
            FN += 1
        else:
            query_subject = query[site]
            truth_subject = truth[site]
            if query_subject[0] != truth_subject[0]:
                FP += 1
            else:
                query_alt = convert(query_subject[1])
                truth_alt = convert(truth_subject[1])
                if query_alt != truth_alt:
                    FP += 1
                else:
                    if set([query_subject[2], query_subject[3]]) != set([truth_subject[2], truth_subject[3]]):
                        FP += 1
                    else:
                        TP += 1


                        if query_subject[4]:
                            TP_phased += 1
                            isdouble, category = get_category(query_subject, min_sv)
                            flag = get_flag(isdouble, category, only_snv, no_sv, no_indel, no_double)
                            if flag:
                                query_record = myrecord(working_chr, site, query_subject[0], set(query_alt), query_subject[5], query_subject[2], query_subject[3], category)
                            
                                if query_subject[5] not in query_blocks:
                                    query_blocks[query_subject[5]] = list()

                                query_blocks[query_subject[5]].append(query_record)
                        
                        if truth_subject[4]:
                            isdouble, category = get_category(truth_subject, min_sv)
                            flag = get_flag(isdouble, category, only_snv, no_sv, no_indel, no_double)
                            if flag:
                                truth_record = myrecord(working_chr, site, truth_subject[0], set(truth_alt), truth_subject[5], truth_subject[2], truth_subject[3], category)
                                if truth_subject[5] not in truth_blocks:
                                    truth_blocks[truth_subject[5]] = list()

                                truth_blocks[truth_subject[5]].append(truth_record)
    
    stat = {"TP": TP, "FP": FP, "FN": FN, "TP_phased": TP_phased}

    return (query_blocks, truth_blocks, stat)








