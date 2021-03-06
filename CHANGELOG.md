Changelog
==========

<!--
Newest changes should be on top.

This document is user facing. Please word the changes in such a way
that users understand how the changes affect the new version.
-->

v0.4-dev
---------------------------
+ Change default blast output from `extend` to `hit`
+ Add option to select type of blast output to use
+ Change hifiasm-flags into a list to ensure parameter value pairs are kept
together
+ Fix a bug where the pipeline crashes if an assembly output is an empty file

v0.3
---------------------------
+ Add support for selecting `hap1`, `hap2` as HiFiasm output
+ Add support for selecting multiple HiFiasm outputs to use
+ Add support for selecting `a_ctg` as HiFiasm output
+ Add option `hifiasm-write-ec`, to set the `--write-ec` option in HiFiasm
+ Include regions from the assembly that overlap, but do not match, the gene of
interest
+ Add separate bam file output for mapped/unmapped contigs
+ Add support for selecting HiFiasm output to use via the `hifiasm-output`
setting
+ Add sample name as prefix to per-gene contigs fasta files

v0.2
---------------------------
+ Reorganise pipeline outputs
+ Output only the best hit to the target genes for each contig
+ Make both 'reference' and 'genes' inputs optional
+ Add functionality for blasting genes of interest against the assembly
+ Add support for multiple bam files, see this [example
configuration](https://github.com/Redmar-van-den-Berg/HiFi-assembly/blob/main/tests/pep/project_config_two_bamfiles.yml)

v0.1
---------------------------
+ Initial release

v0.0.1
---------------------------
+ Initial commit
