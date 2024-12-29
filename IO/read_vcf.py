import logging


def get_idx(header: list, seg: str) -> int:
    default = {"CHROM": 0, "POS": 1, "REF": 3, "ALT": 4, "FORMAT": -2}
    if seg not in header:
        logging.warning(f"{seg} not in vcf header, using default position")
        return default[seg]
    else:
        return header.index(seg)
    


def read_vcf(filename: str, working_chr: str) -> dict:
    # defalut header position
    pos = {"CHROM": 0, "POS": 1, "REF": 3, "ALT": 4, "FORMAT": -2}
    # set variables
    target_chr = "chr" + working_chr
    variants = dict()

    with open(filename, 'r') as vcf_file:
        for line in vcf_file:
            if line.startswith('#') and not(line.startswith("##")):
                header = (line[1: ]).split()
                subject = ["CHROM", "POS", "REF", "ALT", "FORMAT"]
                for i in subject:
                    pos[i] = get_idx(header, i)
            
            if not(line.startswith('#')):
                m = line.split()
                if m[pos["CHROM"]] == target_chr:
                    var_pos = m[pos["POS"]]
                    ref = m[pos["REF"]]
                    alt = m[pos["ALT"]]
                    if (var_pos.isdigit()) and ('.' not in ref) and ('.' not in alt):
                        format_list = m[pos["FORMAT"]].split(':')
                        sample_list = m[-1].split(':')
                        temp_dict = dict(zip(format_list, sample_list))
                        gt_left, gt_right = str(), str()
                        if "GT" in temp_dict:
                            gt = temp_dict["GT"]
                            gt_left, gt_right = gt[0], gt[2]
                            if gt_left != gt_right:
                                isphase = True if gt[1] == '|' else False
                                
                                ps = "Unify"
                                if "PS" in temp_dict:
                                    if temp_dict["PS"] != '.':
                                        ps = temp_dict["PS"]
                                flag = 0
                                if int(var_pos) in variants:
                                    if (variants[int(var_pos)][4] == False):
                                        flag = 1
                                else:
                                    flag = 1

                                if flag:
                                    variants[int(var_pos)] = (ref, alt, gt_left, gt_right, isphase, ps)

                        
    return variants
