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

def get_assembly():
    """ Pick the correct HiFi output file based on the configuration """
    files = list()
    for out in config['hifiasm-output']:
        files.append(rules.assemble.output[out])
    return files

def set_hifiasm_flags():
    """ Set the flags to pass to hifiasm

        This can get a bit complex, since some flags will be set by the
        pipeline, and some can be specified by the user.
    """
    # All the flags specified by the user
    flags = set(config.get('hifiasm-flags', '').split(' '))

    # If the user specified they want error corrected reads, we add the
    # apropriate to the flags for HiFiasm. (This setting also changes the
    # behaviour of the pipeline, which is why this flag isn't passed to HiFiasm
    # in 'hifiasm-flags'
    if 'hifiasm-write-ec' in config:
        flags.add('--write-ec')

    # If we want to use the a_ctg output, we have to pass the --primary flag to
    # HiFiasm
    if 'a_ctg' in config['hifiasm-output']:
        flags.add('--primary')

    # Assign the set of flags back to the configuration
    config['hifiasm-flags'] = flags

class hifiasm():
    """ Namespace for hifiasm functions.

        HiFiasm produces different output files, or the same output files with
        different names, depending on the flags that are passed. To make sure
        snakemake doesn't crash on missing output files, we need helper
        functions that know about the flags for HiFiasm to generate the
        exptected output files
    """
    def r_utg(sample):
        if '--primary' in config['hifiasm-flags']:
            return f'{sample}/assembly/{sample}.r_utg.gfa'
        else:
            return f'{sample}/assembly/{sample}.bp.r_utg.gfa'

    def p_utg(sample):
        if '--primary' in config['hifiasm-flags']:
            return f'{sample}/assembly/{sample}.p_utg.gfa'
        else:
            return f'{sample}/assembly/{sample}.bp.p_utg.gfa'

    def p_ctg(sample):
        if '--primary' in config['hifiasm-flags']:
            return f'{sample}/assembly/{sample}.p_ctg.gfa'
        else:
            return f'{sample}/assembly/{sample}.bp.p_ctg.gfa'

    def a_ctg(sample):
        if '--primary' in config['hifiasm-flags']:
            return f'{sample}/assembly/{sample}.a_ctg.gfa'
        else:
            # Return another file that exists, since multiple 'empty' output
            # files are not supported in Snakemake
            return 'log/{sample}_hifiasm.txt'

    def hap1(sample):
        if '--primary' in config['hifiasm-flags']:
            # Hap1 is not defined when --primary is passed to HiFiasm. Return
            # another file that exists since 'emtpy' output files are not
            # supported in Snakemake
            return '{sample}/assembly/{sample}.ovlp.source.bin'
        else:
            return f'{sample}/assembly/{sample}.bp.hap1.p_ctg.gfa'

    def hap2(sample):
        if '--primary' in config['hifiasm-flags']:
            # Hap1 is not defined when --primary is passed to HiFiasm. Return
            # another file that exists since 'emtpy' output files are not
            # supported in Snakemake
            return '{sample}/assembly/{sample}.ovlp.reverse.bin'
        else:
            return f'{sample}/assembly/{sample}.bp.hap2.p_ctg.gfa'

# Set default values
if 'hifiasm-output' not in config:
    config['hifiasm-output'] = 'r_utg'

# Make sure hifiasm-output is always a list
if isinstance(config['hifiasm-output'], str):
    config['hifiasm-output'] = [config['hifiasm-output']]

# Set hifiasm flags
set_hifiasm_flags()
