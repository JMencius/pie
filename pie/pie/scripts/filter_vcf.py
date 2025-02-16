import os
import sys
import vcf
import re

def is_all_alpha(s):
    return bool(re.fullmatch(r'[A-Za-z]+', s))



def filter_record(record, min_sv : int, no_sex : bool, canonical : bool, only_snv : bool, no_sv : bool, no_indel : bool, no_double : bool, no_centro : bool) -> bool:

    if canonical:
        only_snv = True
        no_double = True
    
    if only_snv:
        for i in record.ALT:
            if len(str(i)) > 1:
                return False
        if len(record.REF) > 1:
            return False


    if no_sv:
        for i in record.ALT:
            if len(str(i)) > min_sv:
                return False
        if len(record.REF) > min_sv:
            return False
    
    if no_indel:
        for i in record.ALT:
            if 1 < len(str(i)) <= min_sv:
                return False
        if 1 < len(record.REF) <= min_sv:
            return False

    if no_double:
        if len(record.ALT) >= 2:
            return False

    if no_centro:
        pass
    
    # some quality control step
    if not(is_all_alpha(str(record.REF))):
        return False

    for i in record.ALT:
        if not(is_all_alpha(str(i))):
            return False

    return True
