
include: 'common.smk'
pepfile: config['pepfile']
config = pep.config.get('HiFi-assembly', dict())
samples = pep.sample_table['sample_name']

rule all:
    input:
        fasta_input = [f'{sample}/{sample}.fasta.gz' for sample in samples],
        assembly = [f'{sample}/{sample}.bp.r_utg.fasta' for sample in samples],
        mapped_contigs = [f'{sample}/{sample}_contigs.bam' for sample in samples],
        descriptions = expand('{sample}/{gene}.tsv', sample=samples, gene=get_genes()),



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

rule map_contigs:
    """ Map the assembled contigs against the reference """
    input:
        contigs = rules.assembly_to_fasta.output,
        reference = config['reference']
    output:
        bam = '{sample}/{sample}_contigs.bam',
        bai = '{sample}/{sample}_contigs.bam.bai'
    log:
        'log/{sample}_contigs.txt'
    container:
        containers['minimap2']
    threads:
        12
    shell: """
        minimap2 -a {input.reference} {input.contigs} -t {threads} 2> {log} \
                | samtools sort -o - >{output.bam}
        samtools index {output.bam}
    """

rule extract_description:
    """ Extract the gene description from the mapped contigs """
    input:
        contigs = rules.map_contigs.output.bam,
        reference = config['reference'],
        script = srcdir('scripts/description-from-bam.py')
    output:
        '{sample}/{gene}.tsv'
    params:
        region = get_region
    log:
        'log/{sample}_{gene}.log'
    container:
        containers['description-extractor']
    shell: """
        python3 {input.script} \
            --bam {input.contigs} \
            --reference {input.reference} \
            --region {params.region} > {output} 2> {log}
    """
