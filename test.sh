#!/bin/bash

python pie.py \
-i ./test-data/HG002/HG002.pacbio.WhatsHap.concated.variants.phased.vcf \
-c ./test-data/HG002/HG002_GRCh38_1_22_v4.2.1_benchmark_hifiasm_v11_phasetransfer.vcf \
-o ./output/HG002 \
--no-sex \
--verbose;
