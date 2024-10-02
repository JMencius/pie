import os
import sys
import vcf

##
#!--------Devlopment log--------!#
## Sep. 1. 2024
## Function to be done read_bed
## fbed filter, no_centro is still in development
##


def read_bed(filename : str) -> dict:
    pass



def filter_record(record, fbed : str, min_sv : int, chrom : set, no_sex : bool, canonical : bool, only_snv : bool, no_sv : bool, no_indel : bool, no_double : bool, no_centro : bool) -> bool:
    # some quality control step
    for i in record.REF:
        if not i.isalpha():
            return False

    for i in record.ALT:
        for j in str(i):
            if not j.isalpha():
                return False

    # read fbed files
    ##PENDING

    if record.CHROM[3:] not in chrom:
        return False

    if no_sex:
        if (record.CHROM[3:] == 'X') or (record.CHROM[3:] == 'Y'):
            return False

    if canonical:
        only_snv = True
        no_double = True
    
    if only_snv:
        for i in record.ALT:
            if len(i) > 1:
                return False
        if len(record.REF) > 1:
            return False


    if no_sv:
        for i in record.ALT:
            if len(i) > min_sv:
                return False
        if len(record.REF) > min_sv:
            return False
    
    if no_indel:
        for i in record.ALT:
            if 1 < len(i) <= min_sv:
                return False
        if 1 < len(record.REF) <= min_sv:
            return False

    if no_double:
        if len(record.ALT) >= 2:
            return False

    if no_centro:
        pass

    return True
