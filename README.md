<img src="https://raw.githubusercontent.com/JMencius/pie/main/pie_logo.png" width="100" alt="PIE logo">

# PIE
PIE (Phasing all-In-one Evaluator) is a command-line tool for evaluating haplotype phasing.

## Installation
Currently, the prebuilt wheel supports Python 3.8 on Linux `x86_64`.
### Option1. Install through pip
Due to name conflict on `PyPI`, we have to use `pie-phasing-eval` for pip installation. It is just a name change, won't affect the command line.
```
conda create -n pie python=3.8;
conda activate pie;
pip install pie-phasing-eval;
```

### Option2. Install through Conda
```
conda create -n pie python=3.8;
conda activate pie;
conda install -c bioconda pie;
```


### Option3. Local installation
```
conda create -n pie python=3.8;
conda activate pie

# Download and extract PIE, then enter the directory containing setup.py
pip install .;
```

## Usage
### Required arguments
| Parameters | Description | Format or example |
|:---:|:---:|:---:| 
| `-i or --input` | file for evaluation | .vcf / .vcf.gz / .bcf |
| `-c or --compare` | ground truth file | .vcf / .vcf.gz / .bcf |
| `-r or --ref` | reference file | .fa / .fasta / .fai |
| `-o or --output` | Output file prefix | exmaple: ./test/output_name |

### Recommend flag
Use `--verbose` to monitor the running process and enable detailed logging. 

### Full usages
```
Usage: pie [OPTIONS]

Options:
  -i, --input TEXT       Input vcf/vcf.gz file for evaluation  [required]
  -n, --name TEXT        User defined sample name, [default: Sample]
  -c, --compare TEXT     Ground truth vcf/vcf.gz file for comparison
                         [required]
  -r, --ref TEXT         Reference file fasta file (.fasta or .fa) or fasta
                         index file (.fai)  [required]
  -o, --output TEXT      Output file prefix,  such as -o ./test/output_name
                         [required]
  -t, --threads INTEGER  Maximum number of parallel threads [default: 24]
  -m, --max-len INTEGER  Maximum variant distance for pairwise calculation
                         [default: 2473538]
  -b, --bed TEXT         .bed file specifying genomic regions to include
                         [default: None]
  --min-sv INTEGER       Minimal length threshold of Structural Variant
                         [default: 30, ALT length > 30 bp is SV]
  --chrom TEXT           Chromosome to evaluate,use comma to join chromosome
                         name e.g. --chrom chr1,chr2,chr3
                         [default:chr1,chr2,chr3,...,chr22]
  --sexchrom TEXT        Sex chromosome,use comma to join chromosome name e.g.
                         --sexchrom chrX,chrY [default: chrX,chrY]
  --mincount INTEGER     Minimum number of phased sites in a phase block
                         [default: 2]
  --block                Output phasing block start and end positions in a BED
                         file
  --no-sex               Ignore sex chromosome
  --canonical            Canonical mode, only evaluate single mutation SNV
                         ignore double heterozygous site
  --only-snv             Only evaluate single nucleotide variation
  --only-indel           Only evaluate insertion and deletion
  --only-sv              Only evaluate structural variant
  --no-snv               Ignore single nucleotide variation
  --no-indel             Ignore insertion and deletion
  --no-sv                Ignore structural variant
  --no-double            Ignore double heterozygous site
  --no-sort              Do not sort chromosome or regions, directly use the
                         input order
  --verbose              Enable verbose mode, printing parameters and progress
                         to standard output
  --version              Show the version and exit.
  --help                 Show this message and exit.
```

## Examples
Suppose `phase.vcf` is the sample VCF file to be evaluated against the ground truth VCF file (`truth.vcf`).
1. (Comprehensive) Evaluate all autosomes
```
pie --verbose -i phase.vcf -c truth.vcf -r ref.fa -o ./output/comprehensive
```

2. (Chromosome specific) Restrict evaluation to selected chromosome
```
pie --verbose --chrom chr6 -i phase.vcf -c truth.vcf -r ref.fa -o ./output/chr6
```

3. (Region specific) Focus only on regions defined in a BED file.
```
pie --verbose -i phase.vcf -c truth.vcf -r ref.fa --bed mhc.bed -o ./output/mhc
```

4. (Filter) Exclude Structural Variants (SV) from the analysis
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

## Test data
Small query and truth test files are provided in [here](https://github.com/JMencius/pie/tree/main/tests) including three formats:
| Filename | Format | Description |
|:---:|:---:|:---:|
| small_query.vcf | VCF | Uncompressed VCF |
| small_query.vcf.gz | VCF.GZ | Compressed VCF |
| small_query.bcf | BCF | Binary VCF |

`small_truth.vcf` is provided as the truth.

The `VCF` format follows the VCF v4.1 specification <https://samtools.github.io/hts-specs/VCFv4.1.pdf>


## Output file
| Output suffix | Description | Condition |
|:---:|:---:|:---:|
| `.variant.stats.csv` | Genotype evaluation result and phased percentage | Always generated |
| `.perchrom.csv` | Per chromosome phasing evaluation result | Always generated |
| `.overall.csv` | Overall sample evaluation result | Always generated |
| `.blocks.bed` | Raw phasing block start and end in `-i` or `--input` file | with `--block` set |

## Resource consumption
`Pie` is expected to complete evaluation within minutes with the default 24 threads on an x86_64 platform.

The actual performance may vary depending on factors such as size of `vcf`, I/O speed, memory speed, and CPU capabilities.


## Acknowledgements
`Pie` is dependent on the following libraries, we are grateful to all the developers/maintainers:
- [click](https://github.com/pallets/click): Python command line
- [cyvcf2](https://github.com/brentp/cyvcf2): VCF/BCF processing
- [pyfastx](https://github.com/lmdu/pyfastx): Reference FASTA processing
- [numba](https://github.com/numba/numba): JIT acceleration
- [sortedcontainers](https://github.com/grantjenks/python-sortedcontainers): Python Sorted Container Types
















