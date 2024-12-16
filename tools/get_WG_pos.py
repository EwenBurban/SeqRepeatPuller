genome_ref=config['genome_ref']# the fasta sequence of the reference genome
TE_bed=config['TE_bed']# the bed file indicating the genomic postion of each TE of interest on the reference genome
TE_refseq=config['TE_refseq']# the fasta sequence of the reference TE
seq_bed=config['seq_bed']# the bed file indicating the position of the sequences of interest on the TE reference sequence
workdir : config['workdir']

rule all:
    input:
        'WG_Xtracted_position.Xbed'

rule get_TE_sequences:
    input:
        genome_ref=genome_ref
        TE_bed=TE_bed
    output:
        temp('TE_OI_sequences.fasta')
    shell:
        """
        bedtools getfasta -fi {input.genome_ref} -bed {input.TE_bed} -fo {output} -s
        """

rule map_TE_sequences_on_TE_refseq:
    input:
        TE_sequences=rules.get_TE_sequences.output,
        TE_refseq=TE_refseq
    output:
        bam=temp('mapped_TE_sequences.bam')
        bam_sorted='mapped_TE_sequences_sorted.bam'
    shell:
        """
        bowtie2-build {input.TE_refseq} bowtie2_index
        ### modify the command to Tessendier parameter
        bowtie2 -x bowtie2_index -U {TE_sequences} | samtools view -bs4 - > {output.bam}
        samtools sort {output.bam} > {output.bam_sorted}
        samtools index {output.bam_sorted}
        """

rule Xtract_WG_position:
    input:
        seq_bed=seq_bed,
        bam_file=rules.map_TE_sequences_on_TE_refseq.output.bam_sorted
    output:
        'WG_Xtracted_position.Xbed'
    shell:
        """
        python3 {binpath}/tools/get_WG_pos_scripts/RefSeqXtractor.py -b {input.bam_file} -p {input.seq_bed} -o {output}
        """






