
include: 'common.smk'
pepfile: config['pepfile']
config = pep.config.get('HiFi-assembly', dict())
samples = pep.sample_table['sample_name']

rule all:
    input:
        fasta_input = [f'{sample}/{sample}.fasta.gz' for sample in samples],
        assembly = [f'{sample}/{sample}.bp.r_utg.fasta' for sample in samples],
        mapped_contigs = [f'{sample}/{sample}_contigs.bam' for sample in samples] if 'reference' in config else [],
        fasta = [f'{sample}/{sample}_contigs_blast.fasta' for sample in samples] if 'genes' in config else [],
        json = [f'{sample}/{sample}_contigs_blast.json' for sample in samples] if 'genes' in config else [],
        xml = [f'{sample}/blastdb/{sample}_genes.xml' for sample in samples] if 'genes' in config else [],

rule bam_to_fasta:
    input:
        get_bamfiles
    output:
        '{sample}/{sample}.fasta.gz'
    log:
        'log/{sample}_bam_to_fasta.txt'
    container:
        containers['minimap2']
    shell: """
        # Make sure the output folder exists
        mkdir -p {wildcards.sample}

        # Remove the output file if it exists, since we will be appending
        rm -f {output} {log}

        for bamfile in {input}; do
            # Gzip files can be appended together
            samtools fasta $bamfile 2>> {log} | gzip >> {output}
        done
    """

rule assemble:
    input:
        fasta = rules.bam_to_fasta.output,
    output:
        # Haplotype-resolved raw unitig graph.
        # This graph keeps all haplotype information.
        r_utg = '{sample}/{sample}.bp.r_utg.gfa',

        # Haplotype-resolved processed unitig graph without small bubbles.
        # Small bubbles might be caused by somatic mutations or noise in data,
        # which are not the real haplotype information.
        p_utg = '{sample}/{sample}.bp.p_utg.gfa'
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
        containers['pyblast']
    shell: """
        python3 {input.script} {input.gfa} {output} 2> {log}
    """

rule map_contigs:
    """ Map the assembled contigs against the reference """
    input:
        contigs = rules.assembly_to_fasta.output,
        reference = config.get('reference', '')
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

rule make_blast_db:
    """ Create a blast database from the assembled contigs """
    input:
        contigs = rules.assembly_to_fasta.output
    params:
        dbname = lambda wildcards, output: output[0][:-4]
    output:
        nhr = '{sample}/blastdb/{sample}.blastdb.nhr',
        nin = '{sample}/blastdb/{sample}.blastdb.nin',
        nsq = '{sample}/blastdb/{sample}.blastdb.nsq'
    log:
        'log/{sample}_make_blast_db.txt'
    container:
        containers['pyblast']
    shell: """
        mkdir -p $(dirname {output.nhr}) 2> {log}

        makeblastdb \
            -input_type fasta \
            -dbtype nucl \
            -in {input} \
            -out {params.dbname} 2>> {log}
    """

rule blast_contigs:
    """ Run blast against the blast database of the contigs """
    input:
        blastdb = rules.make_blast_db.output.nhr,
        genes = config.get('genes', ''),
    params:
        dbname = lambda wildcards, input: input[0][:-4]
    output:
        xml = '{sample}/blastdb/{sample}_genes.xml'
    log:
        'log/{sample}_blast_contigs.txt'
    container:
        containers['pyblast']
    shell: """
        blastn \
            -outfmt 5 \
            -query {input.genes} \
            -db {params.dbname} \
            -out {output} 2> {log}
    """

rule blast_genes:
    """ Blast the specified genes against the assembled contigs """
    input:
        blast_results = rules.blast_contigs.output,
        genes = config.get('genes', ''), # Not actually used in the script
        script = srcdir('scripts/run-blast.py')
    output:
        json = '{sample}/{sample}_contigs_blast.json',
        folder = directory('{sample}/genes'),
        fasta = '{sample}/{sample}_contigs_blast.fasta'
    log:
        'log/{sample}_blast_genes.txt'
    container:
        containers['pyblast']
    shell: """

        # Don't use mkdir -p here, since script appends to output files
        if [ ! -d {output.folder} ]; then
            mkdir {output.folder} 2> {log}
        else
            rm {output.folder}/*.fasta
        fi

        python3 {input.script} \
            --database {input.blast_results} \
            --query {input.genes} \
            --json {output.json} \
            --genes {output.folder} \
            --fasta {output.fasta} 2>> {log}
    """
