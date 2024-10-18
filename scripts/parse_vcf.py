import vcf
import os
import sys
from scripts.filter_vcf import filter_record
from scripts.myrecord import myrecord
from multiprocessing import Pool


def gt_info(gt_tag : str) -> tuple:
    isphased = '|' in gt_tag
    if not isphased:
        return (False, None, False)

    isok = False
    if gt_tag[0].isdigit() and gt_tag[2].isdigit():
        if gt_tag[0] != gt_tag[2]:
            isok = True
    if not isok:
        return (True, None, False)

    isdouble = '2' in gt_tag

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



def process_filter_vcf(filename: str, target_chrom: str, start: int, end: int, min_sv: int, no_sex : bool, canonical : bool, only_snv : bool, no_sv : bool, no_indel : bool, no_double : bool, no_centro : bool) -> dict:

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
    #print(f"read mode is {readmode}")
    
    output = {target_chrom : dict()}

    count = 0
    flag = 0
    for record in vcf_reader:
        # for chromosome parsing
        chr_str = record.CHROM[3:]
        if (chr_str != target_chrom) and (flag == 1):
            break

        if chr_str == target_chrom:
            flag = 1
            for sample in record.samples:
                filter_flag = filter_record(record, min_sv, no_sex, canonical, only_snv, no_sv, no_indel, no_double, no_centro)
                #print(record, filter_flag)
                if filter_flag:
                    # parse vcf only according to GT tag
                    if readmode == "GT":
                        if "GT" in f"{sample.data}":
                            isphased, isdouble, isok = gt_info(sample["GT"])
                        if isok:
                            category = get_category(record.REF, ','.join([str(i) for i in record.ALT]), min_sv)
                            if category != "UNKNOWN":
                                if isphased:
                                    if sample["GT"][0].isdigit and sample["GT"][2].isdigit:
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


                    # parse vcf according to PS tag and GT tag
                    else:
                        if "GT" in f"{sample.data}":
                            isphased, isdouble, isok = gt_info(sample["GT"])
                            if isok:
                                category = get_category(record.REF, ','.join([str(i) for i in record.ALT]), min_sv)
                                if category != "UNKNOWN" and "PS" in f"{sample.data}":
                                        if sample["PS"] != '.' and sample["PS"]:
                                            if isphased:
                                                if sample["GT"][0].isdigit and sample["GT"][2].isdigit:
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


    return output
 

def read_vcf(filename1: str, filename2: str, target: list, min_sv: int, chrom: list, no_sex: bool, canonical: bool, only_snv: bool, no_sv: bool, no_indel: bool, no_double: bool, no_centro: bool, threads: int) -> tuple:
    
    # join target file and target chromosome
    file_chrom = list()
    for f in [filename1, filename2]:
        for c in target:
            file_chrom.append((f, c))
    
    ##print(file_chrom)

    with Pool(threads) as rd:
        temp_read = rd.starmap(process_filter_vcf, [(i, j[0], j[1], j[2], min_sv, no_sex, canonical, only_snv, no_sv, no_indel, no_double, no_centro) for i,j in file_chrom])
    
    file1_assemble = dict()
    file2_assemble = dict()
    for count in range(len(file_chrom)):
        if file_chrom[count][0] == filename1:
            file1_assemble.update(temp_read[count])
        else:
            file2_assemble.update(temp_read[count])


    return (file1_assemble, file2_assemble)



def get_sample_name(vcffile: str) -> str:
    vcf_reader = vcf.Reader(open(vcffile, 'r'))
    return vcf_reader.samples[0]


