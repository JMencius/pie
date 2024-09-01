# pie
Phasing all-In-one Evaluator

## Usage
```
Usage: pie.py [OPTIONS]

Options:
  -i, --input TEXT       Input vcf/bcf file for evaluation  [required]
  -c, --compare TEXT     Truth vcf/bcf file for comparsion
  -r, --ref TEXT         Reference file of the vcf/bcf file
  -o, --output TEXT      Output file  [required]
  -t, --threads INTEGER  Numbers of parallel threads
  -b, --bed TEXT         Bed file to define the centromere region
  --chrom TEXT           Chromosome to evaluate,use comma to connect e.g.
                         --chrom 1,2,3
  --mincount INTEGER     Minimum numbers of phased sites in a phase block
  --no-sex               Ignore sex chromosome
  --only-snv             Only evaluate Single Nucleotide Variation
  --no-sv                Ignore Structural Variant
  --no-indel             Ignore insertion and deletion
  --no-double            Ignore double heterozygous site
  --no-centro            Ignore centromere region
  --only-num             Export only numbers
  --verbose              Verbose mode print some results to stdout
  -v, --version          Print version info and exit
  --help                 Show this message and exit.
```
