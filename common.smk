containers = {
        'python': 'docker://python:latest',
        'hifiasm': 'docker://quay.io/biocontainers/hifiasm:0.16.1--h2e03b76_0',
        # samtools 1.12, minimap 2.20
        'minimap2': 'docker://quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:6203818d4ae215fd8c82ea5f7cc0f50f56066592-0',
        'samtools': 'docker://quay.io/biocontainers/samtools:1.13--h8c37831_0'
}

def get_bamfile(wildcards):
    return pep.sample_table.loc[wildcards.sample, 'bamfile']
