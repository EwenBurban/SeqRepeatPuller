genome_ref=config['genome_ref']# the fasta sequence of the reference genome
te_bed=config['TE_bed']# the bed file indicating the genomic postion of each te of interest on the reference genome
te_refseq=config['TE_refseq']# the fasta sequence of the reference te
seq_bed=config['seq_bed']# the bed file indicating the position of the sequences of interest on the te reference sequence
binpath=config['binpath']
git_commit=config['git_commit']

if 'mode' in config.keys():
    mode=config['mode']
    if mode == "fragmented" :
        fragment_size=config['fragment_size']
        fragment_step=config['fragment_step']
    tag_mode=mode
    print("fragment alignement mode (only fragmented | end-to-end): ", mode)
else :
    mode='end-to-end'
    tag_mode='complete'

if 'clean_xtracted_results' in config.keys():
    cleaning=config['clean_xtracted_results']
    tag_cleaning="_cleaned"
else:
    cleaning='false'
    tag_cleaning=""

rule all:
    input:
        ['wg_xtracted_position.xbed',
        'wg_xtracted_position_metadata.json']

rule get_te_sequences:
    input:
        genome_ref=genome_ref,
        te_bed=te_bed
    threads: 1
    output:
        temp('te_oi_sequences_complete.fasta')
    shell:
        """
        bedtools getfasta -fi {input.genome_ref} -bed {input.te_bed} -fo {output} -s
        """


rule fragment_te_oi_sequences:
    input:
        genome_ref=genome_ref,
        te_bed=te_bed
    threads: 1
    output:
        fasta=temp('te_oi_sequences_fragmented.fasta'),
        frag_bed=temp('fragmented_te_oi.bed')
    shell:
        """
        python3 {binpath}/fragment_TE_bed.py --bed {input.te_bed} \
                --output {output.frag_bed} --fragment_size {fragment_size} --fragment_step {fragment_step}

        bedtools getfasta -fi {input.genome_ref} -bed {output.frag_bed} -fo {output.fasta} -s 
        """



rule map_te_sequences_on_te_refseq:
    input:
        te_sequences=f'te_oi_sequences_{tag_mode}.fasta',
        te_refseq=te_refseq
    threads: 1
    output:
        bam=temp('mapped_te_sequences.bam'),
        bam_sorted=temp('mapped_te_sequences_sorted.bam'),
        bam_index=temp('mapped_te_sequences_sorted.bam.bai')
    shell:
        """
        bowtie2-build {input.te_refseq} bowtie2_index
        bowtie2 -x bowtie2_index -f {input.te_sequences} -N 1 | samtools view -b - > {output.bam}
        samtools sort {output.bam} > {output.bam_sorted}
        samtools index {output.bam_sorted}
        """

rule xtract_wg_position:
    input:
        seq_bed=seq_bed,
        bam_file=rules.map_te_sequences_on_te_refseq.output.bam_sorted,
        bam_index=rules.map_te_sequences_on_te_refseq.output.bam_index
    threads: 1
    output:
        'wg_xtracted_position.xbed'
    shell:
        """
        python3 {binpath}/RefSeqXtractor.py -b {input.bam_file} -p {input.seq_bed} -o {output}
        """

rule clean_xtraction:
    input:
        rules.xtract_wg_position.output
    threads: 1
    output:
       'wg_xtracted_position_cleaned.xbed'
    shell:
        """
        python3 {binpath}/clean_Xtracted_position.py --xbed {input} -m {output}.pkl -o {output}
        """


rule generate_metadata:
    input:
        f'wg_xtracted_position{tag_cleaning}.xbed'
    output:
        'wg_xtracted_position_metadata.json'
    shell:
        """
        python3 {binpath}/generate_metadata.py -i {input} -c "{config}" -o {output}
        """




