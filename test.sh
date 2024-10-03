#!/bin/bash

python pie.py \
-i ./test-data/longphase_HG01109.vcf \
-c ./test-data/HG01109.f1_assembly_v2.dip.vcf \
-o HG01109_longphase_test \
--min-sv 50 \
--only-snv \
--no-double \
--no-sex \
--verbose;
