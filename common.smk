containers = {
        'description-extractor': 'docker://lumc/description_extractor:0.2',
        'python': 'docker://python:latest',
        'hifiasm': 'docker://quay.io/biocontainers/hifiasm:0.16.1--h2e03b76_0',
        # samtools 1.12, minimap 2.20
        'minimap2': 'docker://quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:6203818d4ae215fd8c82ea5f7cc0f50f56066592-0',
        'samtools': 'docker://quay.io/biocontainers/samtools:1.13--h8c37831_0'
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

def get_genes():
    """ Extract the gene names from the configuration """
    # Make sure the order of the genes is always the same
    return sorted(list(config['genes']))

def get_regions():
    """ Extract the gene regions from the configuration.

    In the same order as get_genes()
    """
    return [config['genes'][gene] for gene in get_genes()]

def get_region(wildcards):
    """ Return the region for the gene wildcard """
    return config['genes'][wildcards.gene]
