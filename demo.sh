#!/bin/bash

# Download HG02818 phased VCF
wget -c -O HG02818.nanopore.R10.hac.HapCUT2.concated.variants.phased.vcf.gz 'https://china.scidb.cn/download?fileId=a79b69e3186002aa9651ac63a787bae4'

# Download HG02818 truth VCF
wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG02818/assemblies/year1_f1_assembly_v2_genbank/assembly_qc/dipcall/HG02818.f1_assembly_v2_genbank.dip.vcf.gz

# Download reference file
wget -c https://ont-open-data.s3.amazonaws.com/colo829_2024.03/wf_somatic_variation/sup/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.fai

# Run PIE
pie --verbose -i HG02818.nanopore.R10.hac.HapCUT2.concated.variants.phased.vcf.gz -c HG02818.f1_assembly_v2_genbank.dip.vcf.gz -r GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.fai -o ./HG02818
