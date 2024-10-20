# pie
Phasing all-In-one Evaluator

## Installation
```
conda env create -f pie.yaml;
```

## Usage
```
Usage: pie.py [OPTIONS]

Options:
  -i, --input TEXT       Input vcf/bcf file for evaluation  [required]
  -n, --name TEXT        User defined sample name of the input vcf file,
                         [default: `extract from vcf`]
  -c, --compare TEXT     Truth vcf/bcf file for comparsion  [required]
  -r, --ref TEXT         Reference file [default:GCA_000001405.15_GRCh38_no_al
                         t_analysis_set.fasta]
  -o, --output TEXT      Output tsv file name  [required]
  -t, --threads INTEGER  Maximum numbers of parallel threads
  --fbed TEXT            Bed file to filter out, such as centromere region in
                         data/hg38_centromere.bed
  --min-sv INTEGER       Minimal length of Structral Variant
  --chrom TEXT           Chromosome to evaluate,use comma to connect e.g.
                         --chrom 1,2,3 [default:1-22, X, Y]
  --mincount INTEGER     Minimum numbers of phased sites in a phase block
                         [default: 2]
  --no-sex               Ignore sex chromosome
  --canonical            Canonical mode, only evaluate single mutation SNV
  --only-snv             Only evaluate Single Nucleotide Variation
  --no-sv                Ignore Structural Variant
  --no-indel             Ignore insertion and deletion
  --no-double            Ignore double heterozygous site
  --no-centro            Ignore centromere region
  --verbose              Verbose mode print intermediate results to stdout
  --test                 Run test sample
  --version              Show the version and exit.
  --help                 Show this message and exit.
```
