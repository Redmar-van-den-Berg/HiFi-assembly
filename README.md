[![Continuous Integration](https://github.com/Redmar-van-den-Berg/HiFi-assembly/actions/workflows/ci.yml/badge.svg)](https://github.com/Redmar-van-den-Berg/HiFi-assembly/actions/workflows/ci.yml)
[![PEP compatible](http://pepkit.github.io/img/PEP-compatible-green.svg)](http://pepkit.github.io)
![GitHub release](https://img.shields.io/github/v/release/redmar-van-den-berg/HiFi-assembly)
![Commits since latest release](https://img.shields.io/github/commits-since/redmar-van-den-berg/HiFi-assembly/latest)

# HiFi assembly
Assemble HiFi reads and extract variants for the specified genes.

## Installation
Download the repository from github
```bash
git clone https://github.com/Redmar-van-den-Berg/HiFi-assembly.git
```

Install and activate the
[conda](https://docs.conda.io/en/latest/miniconda.html)
environment.
```bash
conda env create --file environment.yml
conda activate HiFi-assembly
```

## Settings
The settings for this pipeline are defined in a
[PEP](http://pep.databio.org/en/latest/) project file, as shown below.
```yml
pep_version: 2.0.0
sample_table: "samples.csv"
HiFi-assembly:
  hifiasm-flags: "-f0 "
  reference: tests/data/reference/ASL.fasta
  genes:
    ASL: "chr7"
    ASL2: "chr7:1-17758"
```

The samples are defined in a simple csv file.
```csv
sample_name,bamfile
GM24385_ASL,tests/data/GM24385_ASL.bam
```

### Multiple bam files per sample
If you have multiple bam files per sample, you can utilise the
`subsample_table` in the PEP project file, see [this
example](https://github.com/Redmar-van-den-Berg/HiFi-assembly/blob/main/tests/pep/project_config_two_bamfiles.yml)
for details. Be sure to include every sample in the `sample_table`, otherwise
the bam file specified in the `subsample_table` will be ignored.

## How it works
The reads from the bam file(s) are assembled using `HiFiasm` with default
settings, unless you specify additional flags using `hifiasm-flags`. After
assembly, the contigs are mapped against the reference using `minimap2`. Next,
the [mutalyzer description
extractor](https://mutalyzer.nl/description-extractor) is used to directly
determine the difference between the reference genes and the contigs.
Finally, this description is trimmed to exlude insertions and deletions that
occur at the beginning and end of the gene. This way, there are no spurious
differences when the contigs are either larger or smaller than the reference.

## Output files
The following output files are the most relevant
| File path                         | Explanation                             |
| --------------------------------- | --------------------------------------- |
| sample/sample_gene_contigs.bam    | The contigs mapped to the reference     |
| sample/gene.raw.tsv               | The raw descriptions for each contig    |
| sample/gene.trimmed.tsv           | The trimmed descriptions for each contig|
