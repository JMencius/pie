import vcf
import os
import sys
from scripts.filter_vcf import filter_record
from scripts.myrecord import myrecord

##++++++++++++++++++
## Development log
## 2024. Aug. 29
## more quality control step to be added
##++++++++++++++++++





def gt_info(gt_tag : str) -> tuple:
    isphased = '|' in gt_tag
    isdouble = '2' in gt_tag
    if gt_tag[0].isdigit() and gt_tag[2].isdigit():
        isok = True
    else:
        isok = False

    return (isphased, isdouble, isok)



def get_category(ref_tag : str, alt_tag : list, min_sv : int) -> str:
    # quality control of REF and ALT
    # if '*' or other nonsense exist, return UNKNOWN
    flag = 1
    for i in ref_tag + ''.join([str(a) for a in alt_tag]):
        if i.isalpha() or i == ',':
            pass
        else:
            flag = 0
            break
    if flag == 0:
        return "UNKNOWN"
    
    # choose category in SNV, DOUBLE, INDEL, SV
    if len(ref_tag) == 1:
        aflag = 1
        for i in alt_tag:
            if len(i) != 1:
                aflag = 0
        if aflag == 1:
            return "SNV"
    elif len(ref_tag) > min_sv:
        return "SV"
    else:
        for i in alt_tag:
            if len(i) > min_sv:
                return "SV"
    return "INDEL"



def process_filter_vcf(filename: str, fbed: str, min_sv: int, chrom: set, no_sex : bool, canonical : bool, only_snv : bool, no_sv : bool, no_indel : bool, no_double : bool, no_centro : bool) -> tuple:
    chrom_str = {str(i) for i in chrom}
    phase_count = {i : {"SNV": [0, 0], "INDEL": [0, 0], "SV": [0, 0]} for i in chrom_str}
    unphase_count = {i : {"SNV": [0, 0], "INDEL": [0, 0], "SV": [0, 0]} for i in chrom_str}

    output = {i : dict() for i in chrom_str}

    vcf_reader = vcf.Reader(open(filename, 'r'))
    for i in vcf_reader.formats:
        temp = []
        for j in (vcf_reader.formats[i]):
            if j == 'Float':
                temp.append('String')
            else:
                temp.append(j)
        #print(vars(vcf.parser))
        #sys.exit(0)
        vcf_reader.formats[i] = vcf.parser._Format(temp[0], temp[1], temp[2], temp[3])
    
    # determine how to read vcf file according to FORMAT
    readmode = "GT"
    if "GT" not in vcf_reader.formats:
        raise AttributeError("GT tag not in vcf/bcf file format, check the input vcf")
    else:
        if "PS" not in vcf_reader.formats:
            print(f"WARNING PS tag not in vcf/bcf file, the whole chromosome will be treated as a completely phased")
        else:
            readmode = "PS"
    print(f"read mode is {readmode}")


    count = 0
    for record in vcf_reader:
        # for chromosome parsing
        chr_str = record.CHROM[3:]
        if chr_str in chrom_str:
            for sample in record.samples:
                filter_flag = filter_record(record, fbed, min_sv, chrom, no_sex, canonical, only_snv, no_sv, no_indel, no_double, no_centro)
                if filter_flag:
                    # parse vcf only according to GT tag
                    if readmode == "GT":
                        if "GT" in f"{sample.data}":
                            isphased, isdouble, isok = gt_info(sample["GT"])
                        if isok:
                            category = get_category(record.REF, ','.join([str(i) for i in record.ALT]), min_sv)
                            if category != "UNKNOWN":
                                if isphased:
                                    this_record = myrecord(chr_str, 
                                               int(record.POS), 
                                                record.REF, 
                                                set([str(i) for i in record.ALT]), 
                                                "UNIFY", 
                                                sample["GT"][0],
                                                sample["GT"][2],
                                                category)
                                    if "UNIFY" not in output[chr_str]:
                                        output[chr_str]["UNIFY"] = list()
                                    output[chr_str]["UNIFY"].append(this_record)


                                    if isdouble:
                                        phase_count[chr_str][category][1] += 1
                                    else:
                                        phase_count[chr_str][category][0] += 1
                                else:
                                    if isdouble:
                                        unphase_count[chr_str][category][1] += 1
                                    else:
                                        unphase_count[chr_str][category][0] += 1

                
                    # parse vcf according to PS tag and GT tag
                    else:
                        if "GT" in f"{sample.data}":
                            isphased, isdouble, isok = gt_info(sample["GT"])
                            if isok:
                                category = get_category(record.REF, ','.join([str(i) for i in record.ALT]), min_sv)
                                if category != "UNKNOWN":
                                    if isphased:
                                        if isdouble:
                                            phase_count[chr_str][category][1] += 1
                                        else:
                                            phase_count[chr_str][category][0] += 1
                                    else:
                                        if isdouble:
                                            unphase_count[chr_str][category][1] += 1
                                        else:
                                            unphase_count[chr_str][category][0] += 1

                                    if "PS" in f"{sample.data}":
                                        if sample["PS"] != '.' and sample["PS"]:
                                            if isphased:
                                                this_record = myrecord(chr_str,
                                                        int(record.POS),
                                                        record.REF,
                                                        set([str(i) for i in record.ALT]),
                                                        sample["PS"],
                                                        sample["GT"][0],
                                                        sample["GT"][2],
                                                        category)

                                                if sample["PS"] not in output[chr_str]:
                                                    output[chr_str][sample["PS"]] = list()
                                                output[chr_str][sample["PS"]].append(this_record)



    return (output, phase_count, unphase_count)
                        


"""
## test code
if __name__ == "__main__":
    output = process_vcf(sys.argv[1], {1}, 30)
    count = 0
    for i in output['1']["UNIFY"]:
        count += 1
        print(i)
        if count == 10:
            sys.exit(0)
"""
