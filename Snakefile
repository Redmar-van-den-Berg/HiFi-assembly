
include: 'common.smk'
pepfile: config['pepfile']

rule all:
    input:
        fasta_input = expand('{sample}/{sample}.fasta.gz', sample=pep.sample_table['sample_name']),
        assembly = expand('{sample}/{sample}.bp.r_utg.gfa', sample=pep.sample_table['sample_name'])


rule bam_to_fasta:
    input:
        get_bamfile
    output:
        '{sample}/{sample}.fasta.gz'
    log:
        'log/{sample}_bam_to_fasta.txt'
    container:
        containers['samtools']
    shell: """
        mkdir -p {wildcards.sample}
        samtools fasta {input} 2> {log} | gzip > {output}
    """

rule assemble:
    input:
        fasta = rules.bam_to_fasta.output,
    output:
        r_utg = '{sample}/{sample}.bp.r_utg.gfa'
    threads:
        12
    log:
        'log/{sample}_hifiasm.txt'
    container:
        containers['hifiasm']
    shell: """
        hifiasm -o {wildcards.sample}/{wildcards.sample} -t {threads} {input.fasta} > {log}
    """
