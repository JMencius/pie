<img src="pie_logo.png" width = "71" height = "153">

# pie
Phasing all-In-one Evaluator, for haplotype phasing evaluation

## Installation
1. Create new virtual environment
```
conda create -n pie python=3.7;
conda activate pie;
```

2. Navigate to the base directory, which contains `setup.py`. Use `pip` to install `pie`.
```
pip install .;
```


## Usage
### Required argument
| Parameters | Description | Format or example |
|:---:|:---:|:---:| 
| `-i or --input` | file for evaluation | .vcf / .vcf.gz |
| `-c or --compare` | ground truth file | .vcf / .vcf.gz |
| `-r or --ref` | reference file | .fa / .fasta / .fai |
| `-o or --output` | Output file prefix | exmaple: ./test/output_name |

### Recommend flag
Use `--verbose` to monitor the running process and enable detailed logging.

```
Usage: pie [OPTIONS] 

Options:
Input and output files argument:
  -i, --input TEXT       Input vcf/vcf.gz file for evaluation  [required]
  -n, --name TEXT        User defined sample name, [default: Sample]
  -c, --compare TEXT     Ground truth vcf/vcf.gz file for comparison  [required]
  -r, --ref TEXT         Reference file fasta file (.fasta or .fa) or fasta index file (.fai)  [required]
  -o, --output TEXT      Output file prefix,  such as -o ./test/output_name  [required]
  -t, --threads INTEGER  Maximum numbers of parallel threads [default: 24]

Pairwise length argument:
  -m, --max-len INTEGER  Maximum variant distance for pairwise calculation [default: 250000]

Specific site or chromosome, SV definiton:
  -b, --bed TEXT         .bed file specifying genomic regions to include.
  --min-sv INTEGER       Minimal length threshold of Structral Variant [default: 30], > 30 bp is SV
  --chrom TEXT           Chromosome to evaluate,use comma to join chromosome name e.g. --chrom chr1,chr2,chr3 [default:chr1,chr2,chr3,...,chr22]
  --sexchrom TEXT        Sex chromosme,use comma to join chromosome name e.g. --sexchrom chrX,chrY [default: chrX,chrY]
  --mincount INTEGER     Minimum numbers of phased sites in a phase block [default: 2]

flag:
  --block                Output phasing block start and end positions in a BED file
  --no-sex               Ignore sex chromosome
  --canonical            Canonical mode, only evaluate single mutation SNV ignore double heterozygous site
  --only-snv             Only evaluate single nucleotide variation
  --only-indel           Only evaluate insertion and deletion
  --only-sv              Only evaluate structural variant
  --no-snv               Ignore single nucleotide variation
  --no-indel             Ignore insertion and deletion
  --no-sv                Ignore structural variant
  --no-double            Ignore double heterozygous site
  --verbose              Enable verbose mode, printing parameters and progress to standard output
  --version              Show the version and exit.
  -h, --help             Show this message and exit.
```

## Exmaples
Supoose `phase.vcf` is the sample VCF file to be evaluated against the ground truth VCF file (`truth.vcf`).
1. (Comprehensive) Evaluate all autosome
```
pie --verbose -i phase.vcf -c truth.vcf -r ref.fa -o ./output/comprehensive
```

2. (Chromosome specific) Restrict evaluation to selected chromosome
```
pie --verbose -i phase.vcf -c truth.vcf -r ref.fa -o ./output/chr6
```

3. (Region specific) Focus only on regions defined in a BED file.
```
pie --verbose -i phase.vcf -c truth.vcf -r ref.fa --bed mhc.bed -o ./output/mhc
```

4. (Filter) Exclude structural variants (SV) from the analysis
```
pie --verbose -i phase.vcf -c truth.vcf -r ref.fa --no-sv -o ./output/no_sv
```

5. (More output) Output raw phasing block start and end positions to a BED file
```
pie --verbose -i phase.vcf -c truth.vcf -r ref.fa --block -o ./output/more
```

6. (Miscellaneous) Evaluate only SNVs on chromosome 6 and output raw block start and end positions to a BED file.
```
pie --verbose -i phase.vcf -c truth.vcf -r ref.fa --chrom chr6 --only-snv --block -o ./output/misc
```

## Output file
| Output suffix | Description | Condition |
|:---:|:---:|:--:|
| `.variant.stats.csv` | Genotype evaluation result and phased percentage | Always generated |
| `.perchrom.csv` | Per chromosome phasing evaluation result | Always generated |
| `.overall.csv` | Overall sample evaluation result | Always generated |
| `.blocks.bed` | Raw phasing block start end in `-i` or `--input` file | with `--block` set |


## Resouce consumption
`Pie` is expected to completed evaluation within minutes.
The actual performance may vary depending on factors such as size of `vcf`, I/O speed, memory speed, and CPU capabilities.














