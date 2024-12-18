

import os
import re
import sys
#### mandatory arguments
input_=config['input_']
format_file= config['format_file']
pos_file=config['pos_file']
binpath=config['binpath']
#### optional arguments
if 'genome_folder' in config.keys():
    genome_folder=config['genome_folder']

if 'read_naming' in config.keys() :
    readTrack_naming_convention=config['readTrack_naming']
else : 
    readTrack_naming_convention='_R1|_R2'
R1_pattern,R2_pattern=readTrack_naming_convention.split('|')

if 'aligner' in config.keys():
    aligner=config['aligner']
else :
    aligner='bowtie2'

if 'fastp_param' in config.keys() :
    fastp_param=config['fastp_param']
else:
    fastp_param=''

## retrieve name tag list and input path

if os.path.exists(input_):  # check for input_ exists, then it'll check if 
    # input_ is a file or a directory. If it's a directory, it'll take all file
    # that end with the file format provided. Otherwise, it'll take only the file 
    # input_ as input.
    if os.path.isfile(input_):  
        print("It is a normal file.")
        print('input_ detected as a single file')
        sample_list = input_
        print(f'input file considered: {input_}')
        sample_name = re.sub('.*/','',sample_list)
        name_tags = sample_name.replace('.' + format_file,'')
        input_path = input_.replace(sample_list,'')
    elif os.path.isdir(input_): 
        print("It is a directory.")
        print('input_ detected as a directory, file format provided will be considered')
        sample_list = [os.path.join(input_,f) for f in os.listdir(input_) if f.endswith(format_file)]
        sample_list_joined = '\n'.join(sample_list)
        print(f'input files considered: \n{sample_list_joined}')
        sample_name = [f for f in os.listdir(input_) if f.endswith(format_file)]
        name_tags = [f.replace('.' + format_file,'') for f in sample_name]
        input_path=input_
    else:
        print("Something went wrong !! The provided file is a special file (socket, FIFO, device file), which is not covered by SeqRepeatPuller")
        sys.exit(0)
else:
    print("The input path provided do not exist")
    sys.exit(0)

## Determine from format file the output list
output_files=['Xtracted_seq.txt','sorted.bam','sorted.bam.bai']
dir1=input_path
dir2=input_path
dir_final=workdir
if format_file!='bam' :
    if aligner=='bismark':
        output_files.append('bismark.bam')
    elif aligner=='star':
        output_files.append('star.bam')
    elif aligner=='bowtie2':
        output_files.append('bowtie2.bam')
    

    if format_file!='merged.fastq.gz':
        output_files.append('merged.fastq.gz')
        name_tags=list(set([re.sub(read_naming_convention,'',x) for x in name_tags]))
            
    readTrack_names=[R1_pattern,R2_pattern]
    dir1 = input_path
    dir2 = workdir
    dir_final = workdir
else:
    readTrack_names=''

if format_file not in ['fq', 'fastq', 'bam']:
    print("format file provided non allowed. Only bam or fastq and fq are allowed")
    sys.exit(0)


wildcard_constraints:
    aligner=aligner,
    name_tag='|'.join(name_tags),
    dir1=dir1,
    dir2=dir2,
    dir_final=dir_final

rule all:
    input:
 #       inputs=expand('{input_path}/{name_tag}{readtrack_name}.{format_file}',input_path=input_path,readtrack_name=readtrack_names,name_tag=name_tag,format_file=format_file),
        outputs=expand('{name_tags}_{output_files}',name_tags=name_tags,output_files=output_files)
    shell:
        """
        """


rule merge:
    input:
        r1='{dir1}/{name_tag}{r1_pattern}.{format_file}',
        r2='{dir1}/{name_tag}{r2_pattern}.{format_file}'
    params:
        fastp_param
    output:
        '{dir_final}/{name_tag}_merged.fastq.gz'
    shell:
        """
        a="{format_file}"
        fastp -i {input.r1} -m -i {input.r2} {params} --merged_out {output}

        """

rule mapping_bismark:
    input:
        '{name_tag}_merged.fastq.gz'
    output:
        bam='{name_tag}_bismark.bam',
        tmp_dir=temp(directory('{name_tag}_bismark_tempdir')),
        tmp_uniq_bam=temp('{name_tag}_bismark_tmp.bam'),
        tmp_random_bam=temp('{name_tag}_bismark_tmp.ambig.bam')
    shell:
        """
        mkdir -p {output.tmp_dir}/
        bismark \
        --fasta \
        --ambig_bam \
        --ambiguous \
        --unmapped \
        --gzip \
        --rg_tag \
        --rg_id {wildcards.name_tag} \
        --rg_sample {wildcards.name_tag} \
        --non_bs_mm \
        --parallel 1 \
        --temp_dir {output.tmp_dir} \
        --basename {wildcards.name_tag}_bismark_tmp \
        -n 1 \
        -l 15 \
        --genome_folder {genome_folder} \
        --single_end {input}
        samtools merge -o {output.bam} {output.tmp_uniq_bam} {output.tmp_random_bam} 
        """

rule mapping_star:
    input:
        '{name_tag}_merged.fa'
    output:
        '{name_tag}_star.bam'
    shell:
        """
        """

rule mapping_bowtie2:
    input:
        '{name_tag}_merged.fa'
    output:
        '{name_tag}_bowtie2.bam'
    shell:
        """

        """

rule mapping:
    input:
        '{dir_final}/{name_tag}_merged.fa'
    output:
        '{dir_final}/{name_tag}.bam'
    run:
        if aligner == 'bowtie2':
            shell("
            bowtie2 -x {genome_folder} -U {input} | samtools view -bs4 - > {output}
                  ")

        elif aligner == 'bismark':
            shell("
        mkdir -p {dir_final}/.tmp_{wildcards.name_tag}/
        bismark \
        --fasta \
        --ambig_bam \
        --ambiguous \
        --unmapped \
        --gzip \
        --rg_tag \
        --rg_id {wildcards.name_tag} \
        --rg_sample {wildcards.name_tag} \
        --non_bs_mm \
        --parallel 1 \
        --temp_dir {dir_final}/.tmp_{wildcards.name_tag}/ \
        --basename {dir_final}/{wildcards.name_tag}_bismark \
        -n 1 \
        -l 15 \
        --genome_folder {genome_folder} \
        --single_end {input}
        samtools merge -o {output} {dir_final}/{wildcards.name_tag}_bismark.bam {dir_final}/{wildcards.name_tag}_bismark.ambig.bam")

        elif aligner == 'start':
            shell('touch {output}')

rule sort_and_index:
    input:
        '{dir2}/{name_tag}.bam'
    output:
        sortbam=temp('{dir_final}/{name_tag}_sorted.bam'),
        index=temp('{dir_final}/{name_tag}_sorted.bam.bai')
    shell:
        """
        samtools sort -o {output.sortbam} {input}
        samtools index {output.sortbam}
        """


rule xtraction:
    input:
        bam='{dir_final}/{name_tag}_sorted.bam',
        index='{dir_final}/{name_tag}_sorted.bam.bai',
        pos_file=pos_file
    output:
        '{dir_final}/{name_tag}_Xtracted_seq.txt'
    shell:
        """
        python3 {binpath}/core/xtractor.py -b {input.bam} -p {input.pos_file} -o {output}  
        """

rule generate_metadata:
    input:
        '{name_tag}_xtracted_seq.txt'
    output:
        temp('{name_tag}_metadata.json')
    run:
        import json



        
rule add_metadata_and_zip:
    input:
        xtract='{name_tag}_xtracted_seq.txt',
        metadata='{name_tag}_metadata.json'
    output:
        'processeddata/{name_tag}_xtracted_seq.zip'
    shell:
        """
        zip {output} {input.xtract} {input.metadata}
        """








