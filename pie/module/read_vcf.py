import logging
from cyvcf2 import VCF
from pie_class import variant


def get_type(ref: str, alt: list, min_sv: int) -> tuple:
    if len(alt) > 1:
        double = True
    else:
        double = False

    max_len = max(len(ref), max([len(i) for i in alt]))
    if max_len == 1:
        vtype = "SNV"
    elif 1 < max_len <= min_sv:
        vtype = "INDEL"
    else:
        vtype = "SV"

    return (double, vtype)



def read_vcf(filename: str, working_chr: str, bed_target: str, min_sv: int) -> dict:
    chr_variants = dict()

    for variant in VCF(filename, threads = 1):
        if variant.CHROM == working_chr:
            isphased = varaint.gt_phases[0]
            ps_tag = str(variant.format('PS')[0])
            left, right, _ = variant.genotypes[0]
            if left != right:
                isdouble, variant_type = get_type(variant.REF, variant.ALT, min_sv)
                current_variant = variant(variant.CHROM, variant.POS, variant.REF, variant.ALT, ps_tag, left, right, variant_type, isphased, isdouble)
            
                chr_variants[variant.POS] = current_variant
    
    if not bed_target:
        return chr_variants
    else:
        regions = bed_target[working_chr]
        region_variants_dict = dict()
        for i in regions:
            this_region = dict()
            for site in range(i[0], i[1] + 1):
                if site in chr_variants:
                    this_region[site] = chr_variants[site]
            region_variants_dict[i] = this_region

        # free memory
        del chr_variants

        return region_variants_dict


