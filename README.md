# pie
Phasing all-In-one Evaluator

## Installation
1. Create new virtual environment
```
conda create -n pie python=3.7;
conda activate pie;
```

2. Navigate to the `README.md` directory
```
pip install .;
```


## Usage
```
Usage: pie [OPTIONS]

Options:
  -i, --input TEXT       Input vcf/bcf file for evaluation  [required]
  -n, --name TEXT        User defined sample name of the input vcf file,
                         [default: Sample]
  -c, --compare TEXT     Truth vcf/bcf file for comparsion  [required]
  -r, --ref TEXT         Reference file fasta or fasta.fai  [required]
  -o, --output TEXT      Output tsv file prefix, path can be added before the
                         prefix, such as -o /test/output_name  [required]
  -t, --threads INTEGER  Maximum numbers of parallel threads
  --bed TEXT             Regions to only include, defined in bed file
  --block TEXT           Output phasing block start and end to bed file, such
                         as --block blocks.bed
  --min-sv INTEGER       Minimal length of Structral Variant
  --chrom TEXT           Chromosome to evaluate,use comma to join chromosome
                         name e.g. --chrom chr1,chr2,chr3
                         [default:chr1-chr22,]
  --sexchrom TEXT        Sex chromosme,use comma to join chromosome name e.g.
                         --sexchrom chrX,chrY [default:chrX,chrY,]
  --mincount INTEGER     Minimum numbers of phased sites in a phase block
                         [default: 2]
  --no-sex               Ignore sex chromosome
  --canonical            Canonical mode, only evaluate single mutation SNV
  --only-snv             Only evaluate Single Nucleotide Variation
  --no-sv                Ignore Structural Variant
  --no-indel             Ignore insertion and deletion
  --no-double            Ignore double heterozygous site
  --verbose              Verbose mode print intermediate results to stdout
  --version              Show the version and exit.
  --help                 Show this message and exit.
```
