import logging
from cyvcf2 import VCF
from pie.module.pie_class import pievariant
import sys



def get_type(ref: str, alt: list, min_sv: int, left: int, right: int) -> tuple:
    if len(alt) > 1:
        double = True
    else:
        double = False
    
    variant_pool = list(ref) + alt
    if not ((0 <= left < len(variant_pool)) and (0 <= right < len(variant_pool))): 
        return None

    max_len = max(len(ref), len(variant_pool[left]), len(variant_pool[right]))
    if max_len == 1:
        vtype = "SNV"
    elif 1 < max_len <= min_sv:
        vtype = "INDEL"
    else:
        vtype = "SV"

    return (double, vtype)



def read_vcf(filename: str, working_chr: str, bed_target: dict, min_sv: int) -> dict:
    chr_variants = dict()
    for variant in VCF(filename, threads = 1):
        if variant.CHROM == working_chr:
            isphased = variant.gt_phases[0]
            if "PS" in variant.FORMAT:
                if variant.format("PS"):
                    ps_tag = str(variant.format("PS")[0])
                else:
                    ps_tag = "UNK"
            else:
                ps_tag = "UNK"
            left, right, _ = variant.genotypes[0]
            if (left != right) and ('.' not in variant.gt_bases[0]):
                type_result = get_type(variant.REF, variant.ALT, min_sv, left, right)
                if not type_result:
                    logging.warning(f"GT tag exceed length in {working_chr} position {variant.POS}, which will be ignored")
                    continue
                
                isdouble, variant_type = type_result
                
                current_variant = pievariant(variant.CHROM, variant.POS, variant.REF, variant.ALT, ps_tag, left, right, variant_type, isphased, isdouble)
                
                flag = True
                if variant.POS in chr_variants:
                    if chr_variants[variant.POS].isphased:
                        flag = False
                if flag:
                    chr_variants[variant.POS] = current_variant
    if not bed_target:
        return chr_variants
    else:
        region_variants_dict = dict()
        if working_chr in bed_target:
            regions = bed_target[working_chr]
            for i in regions:
                this_region = dict()
                for site in range(i[0], i[1] + 1):
                    if site in chr_variants:
                        this_region[site] = chr_variants[site]
                region_variants_dict[(working_chr, i[0], i[1])] = this_region

            # free memory
            del chr_variants

        return region_variants_dict


