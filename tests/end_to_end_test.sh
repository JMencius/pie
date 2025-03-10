pie \
-t 1 \
-i ./small_query.vcf \
-c ./small_truth.vcf \
--chrom chr1 \
-r ./GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.fai \
-o ./test \
--verbose;
