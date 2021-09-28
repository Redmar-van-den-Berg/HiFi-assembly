containers = {
        'hifiasm': 'docker://quay.io/biocontainers/hifiasm:0.16.1--h2e03b76_0',
        'samtools': 'docker://quay.io/biocontainers/samtools:1.13--h8c37831_0'
}

def get_bamfile(wildcards):
    return pep.sample_table.loc[wildcards.sample, 'bamfile']
