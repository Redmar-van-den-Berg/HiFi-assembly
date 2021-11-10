pepfile: config['pepfile']
config = pep.config.get('HiFi-assembly', dict())
samples = pep.sample_table['sample_name']
include: 'common.smk'

rule all:
    input:
        assembly = [f'{sample}/assembly/{sample}.bp.{config["hifiasm-output"]}.fasta' for sample in samples],
        mapped_contigs = [f'{sample}/bamfile/{sample}_contigs.bam' for sample in samples] if 'reference' in config else [],
        json = [f'{sample}/blast/{sample}_contigs_blast.json' for sample in samples] if 'genes' in config else [],
        ec_bam = [f'{sample}/bamfile/{sample}.ec.bam' for sample in samples] if 'hifiasm-write-ec' in config and 'reference' in config else [],

rule bam_to_fasta:
    input:
        get_bamfiles
    output:
        '{sample}/input/{sample}.fasta.gz'
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
        r_utg = '{sample}/assembly/{sample}.bp.r_utg.gfa',

        # Haplotype-resolved processed unitig graph without small bubbles.
        # Small bubbles might be caused by somatic mutations or noise in data,
        # which are not the real haplotype information.
        p_utg = '{sample}/assembly/{sample}.bp.p_utg.gfa',

        # Assembly graph of primary contigs. This graph includes a complete
        # assembly with long stretches of phased blocks.
        # NOTE: The assembly contains occasional phase switching, use with
        # extreme caution
        p_ctg = '{sample}/assembly/{sample}.bp.p_ctg.gfa',

        # Alternate assembly contig graph (prefix.a_ctg.gfa). This graph
        # consists of all assemblies that are discarded in primary contig graph
        a_ctg = '{sample}/assembly/{sample}.bp.a_ctg.gfa',

        # Error corrected input reads
        ec_fasta = '{sample}/assembly/{sample}.ec.fa' if 'hifiasm-write-ec' in config else [],
    params:
        flags = config.get('hifiasm-flags', ''),
        write_ec = '--write-ec' if 'hifiasm-write-ec' in config else '',
    threads:
        12
    log:
        'log/{sample}_hifiasm.txt'
    container:
        containers['hifiasm']
    shell: """
        hifiasm -o {wildcards.sample}/assembly/{wildcards.sample} \
        -t {threads} {params.flags} {params.write_ec} \
        {input.fasta} 2> {log}
    """

rule assembly_to_fasta:
    """ Convert the gfa assembly graph to fasta format """
    input:
        gfa = get_assembly(),
        script = srcdir('scripts/gfa-to-fasta.py')
    output:
        f'{{sample}}/assembly/{{sample}}.bp.{config["hifiasm-output"]}.fasta'
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
        mapped_bam = '{sample}/bamfile/{sample}_contigs.bam',
        mapped_bai = '{sample}/bamfile/{sample}_contigs.bam.bai',
        unmapped_bam = '{sample}/bamfile/{sample}_contigs_unmapped.bam',
        unmapped_bai = '{sample}/bamfile/{sample}_contigs_unmapped.bam.bai',
    log:
        minimap = 'log/{sample}_contigs_minimap.txt',
        samtools = 'log/{sample}_contigs_samtools.txt',
    container:
        containers['minimap2']
    threads:
        12
    shell: """
        minimap2 -a {input.reference} {input.contigs} -t {threads} 2> {log.minimap} \
                | samtools sort \
                | samtools view \
                    -h \
                    -b \
                    -o {output.mapped_bam} \
                    -U {output.unmapped_bam} \
                    -F 4 2>{log.samtools}
                samtools index {output.mapped_bam}
                samtools index {output.unmapped_bam}
    """

rule map_ec_reads:
    """ Map the assembled contigs against the reference """
    input:
        contigs = rules.assemble.output.ec_fasta,
        reference = config.get('reference', '')
    output:
        mapped_bam = '{sample}/bamfile/{sample}.ec.bam',
        mapped_bai = '{sample}/bamfile/{sample}.ec.bam.bai',
        unmapped_bam = '{sample}/bamfile/{sample}.ec.unmapped.bam',
        unmapped_bai = '{sample}/bamfile/{sample}.ec.unmapped.bam.bai',
    log:
        minimap = 'log/{sample}_map_ec_reads.txt',
        samtools = 'log/{sample}_map_ec_reads_samtools.txt',
    container:
        containers['minimap2']
    threads:
        12
    shell: """
        minimap2 -a {input.reference} {input.contigs} -t {threads} 2> {log.minimap} \
                | samtools sort \
                | samtools view \
                    -h \
                    -b \
                    -o {output.mapped_bam} \
                    -U {output.unmapped_bam} \
                    -F 4 2>{log.samtools}
                samtools index {output.mapped_bam}
                samtools index {output.unmapped_bam}
    """

rule make_blast_db:
    """ Create a blast database from the assembled contigs """
    input:
        contigs = rules.assembly_to_fasta.output,
    params:
        dbname = lambda wildcards, output: output[0][:-4]
    output:
        nhr = '{sample}/blast/{sample}.blastdb.nhr',
        nin = '{sample}/blast/{sample}.blastdb.nin',
        nsq = '{sample}/blast/{sample}.blastdb.nsq',
        folder = directory('{sample}/blast'),
    log:
        'log/{sample}_make_blast_db.txt'
    container:
        containers['pyblast']
    shell: """
        mkdir -p {output.folder} 2> {log}

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
        xml = '{sample}/blast/{sample}_blast.xml'
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

rule parse_blast_results:
    """ Parse the blast results of the genes against the assembled contigs """
    input:
        blast_results = rules.blast_contigs.output,
        contigs = rules.assembly_to_fasta.output,
        genes = config.get('genes', ''), # Not actually used in the script
        folder = rules.make_blast_db.output.folder,
        script = srcdir('scripts/parse-blast.py')
    output:
        json = '{sample}/blast/{sample}_contigs_blast.json',
    log:
        'log/{sample}_parse_blast_result.txt'
    container:
        containers['pyblast']
    shell: """

        # Make sure there are no fasta files, since script appends to output files
        rm -f {input.folder}/*.fasta 2> {log}

        python3 {input.script} \
            --database {input.blast_results} \
            --query {input.genes} \
            --json {output.json} \
            --genes {input.folder} \
            --gene-prefix {wildcards.sample} \
            --contigs {input.contigs} \
            2>> {log}
    """
