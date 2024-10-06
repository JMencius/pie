#!/bin/bash

python pie.py \
-i ./test-data/longphase_HG01109.vcf \
-c ./test-data/HG01109.f1_assembly_v2.dip.vcf \
-o ~/test/HG01109_longphase_test \
--canonical \
--no-sex \
--verbose;
