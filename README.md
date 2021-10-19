[![Continuous Integration](https://github.com/Redmar-van-den-Berg/HiFi-assembly/actions/workflows/ci.yml/badge.svg)](https://github.com/Redmar-van-den-Berg/HiFi-assembly/actions/workflows/ci.yml)
[![PEP compatible](http://pepkit.github.io/img/PEP-compatible-green.svg)](http://pepkit.github.io)
![GitHub release](https://img.shields.io/github/v/release/redmar-van-den-berg/HiFi-assembly)
![Commits since latest release](https://img.shields.io/github/commits-since/redmar-van-den-berg/HiFi-assembly/latest)

# HiFi assembly
Assemble HiFi reads and extract the specified genes from the assembly.

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
  # This hifias-flag enables the 'low memory' mode, usefull for testing
  hifiasm-flags: "-f0 "
  reference: tests/data/reference/ASL.fasta
  genes: tests/data/reference/ASL.fasta
```

The samples are defined in a simple csv file.
```csv
sample_name,bamfile
GM24385,tests/data/GM24385_ASL.bam
```
### Supported settings
The following settings are available for the pipeline, place them under the
`HiFi-assembly` section in the project configuration.
| Option                            | Type              | Explanation                             |
| --------------------------------- | ----------------- | --------------------------------------- |
| reference                         | Optional file     | If specified, the contigs will be mapped to the reference |
| genes                             | Optional file     | If specified, the genes will be compared to the contigs using BLAST |
| hifiasm-flags                     | Optional string   | Flags to pass to HiFiasm                |

### Multiple bam files per sample
If you have multiple bam files per sample, you can utilise the
`subsample_table` in the PEP project file, see [this
example](https://github.com/Redmar-van-den-Berg/HiFi-assembly/blob/main/tests/pep/project_config_two_bamfiles.yml)
for details. Be sure to include every sample in the `sample_table`, otherwise
the bam file specified in the `subsample_table` will be ignored.

## How it works
### Assembly
The reads from the bam file(s) are assembled using `HiFiasm` with default
settings. You can control the behaviour of the assembly using the
`hifiasm-flags` in the project configuration file.
This pipeline uses the haplotype-resolved raw unitig graph from HiFiasm, since
this graph [contains all haplotype
information](https://hifiasm.readthedocs.io/en/latest/interpreting-output.html).

The assembly is placed in the `sample/assembly` folder.

### Align contigs to the reference
If a `reference` has been specified in the project configuration file, the
contigs are mapped to the reference using `minimap2`. This allows for visual
inspection of the contigs in IGV. Additionally, unexpected results such as
unmapped contigs, or contigs with large scale deletions can be identified from
the bam file.

The bam file is place in the `sample/bamfile` folder.

### Blast genes of interest against the contigs
If a `genes` FASTA file has been specified, these will be blasted against the
contigs to assign them to the corresponding genes.

The full blast results in XML format are placed in
`sample/blast/sample_blast.xml`. Additionally, the section that matches the
sequence in the `genes` FASTA file will be placed in a FASTA file with the
corresponding gene name. For example, if contig `utg000010l` contains the
`CYP2D6` gene, the sequence content of `utg000010l` that matches `CYP2D6` will
be placed in `sample/blast/CYP2D6.fasta`. The fasta header will include
information about the region of the contig that matches, i.e.
`utg000010l:9807-6104 (CYP2D6)`.
