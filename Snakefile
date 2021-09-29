
include: 'common.smk'
pepfile: config['pepfile']
config = pep.config.get('HiFi-assembly', dict())

rule all:
    input:
        fasta_input = expand('{sample}/{sample}.fasta.gz', sample=pep.sample_table['sample_name']),
        assembly = expand('{sample}/{sample}.bp.r_utg.fasta', sample=pep.sample_table['sample_name'])


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
    params:
        flags = config.get('hifiasm-flags', '')
    threads:
        12
    log:
        'log/{sample}_hifiasm.txt'
    container:
        containers['hifiasm']
    shell: """
        hifiasm -o {wildcards.sample}/{wildcards.sample} \
        -t {threads} {params} \
        {input.fasta} 2> {log}
    """

rule assembly_to_fasta:
    """ Convert the gfa assembly graph to fasta format """
    input:
        gfa = rules.assemble.output.r_utg,
        script = srcdir('scripts/gfa-to-fasta.py')
    output:
        '{sample}/{sample}.bp.r_utg.fasta'
    log:
        'log/{sample}_assembly_to_fasta.txt'
    container:
        containers['python']
    shell: """
        python3 {input.script} {input.gfa} {output} 2> {log}
    """
