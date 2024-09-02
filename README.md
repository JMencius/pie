# pie
Phasing all-In-one Evaluator

## Usage
```
Usage: pie.py [OPTIONS]

Options:
  -i, --input TEXT       Input vcf/bcf file for evaluation  [required]
  -c, --compare TEXT     Truth vcf/bcf file for comparsion  [required]
  -r, --ref TEXT         Reference file [default:
                         GCA_000001405.15_GRCh38_no_alt]
  -o, --output TEXT      Output [default: Stdout]
  -t, --threads INTEGER  Maximum numbers of parallel threads
  --fbed TEXT            Bed file to filter out, such as centromere region in
                         data/hg38_centromere.bed
  --min-sv INTEGER       Minimal length of Structral Variant
  --chrom TEXT           Chromosome to evaluate,use comma to connect e.g.
                         --chrom 1,2,3
  --mincount INTEGER     Minimum numbers of phased sites in a phase block
  --no-sex               Ignore sex chromosome
  --canonical            Canonical mode, only evaluate single mutation SNV
  --only-snv             Only evaluate Single Nucleotide Variation
  --no-sv                Ignore Structural Variant
  --no-indel             Ignore insertion and deletion
  --no-double            Ignore double heterozygous site
  --no-centro            Ignore centromere region
  --only-num             Export only numbers
  --verbose              Verbose mode print some results to stdout
  --version              Show the version and exit.
  --help                 Show this message and exit.
```
