containers = {
        # biopython 1.78, blast 2.9.0, pysam
        'pyblast': 'docker://quay.io/biocontainers/mulled-v2-f55f1045089e924c317dd7bb1e0cf2a139bdb66e:2f42a30ae9e34e2fcbe8ef235ab723b0a279b555-0',
        'hifiasm': 'docker://quay.io/biocontainers/hifiasm:0.16.1--h2e03b76_0',
        # samtools 1.12, minimap 2.20
        'minimap2': 'docker://quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:6203818d4ae215fd8c82ea5f7cc0f50f56066592-0',
}

def get_bamfiles(wildcards):
    """ Return a list of bamfiles for sample """
    bamfiles = pep.sample_table.loc[wildcards.sample, 'bamfile']

    # If a single bamfile is specified, bamfiles will be a string
    if isinstance(bamfiles, str):
        return [bamfiles]
    # If multiple bam files were specified, bamfiles will be a list
    else:
        return bamfiles
