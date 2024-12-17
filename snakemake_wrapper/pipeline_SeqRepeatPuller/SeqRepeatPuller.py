

import os
import re
input_=config['input_']
format_file=config['format_file']
genome_folder=config['genome_folder']

if 'read_naming' in config.keys() :
    readTrack_naming_convention=config['readTrack_naming']
else : 
    readTrack_naming_convention='_R1|_R2'
R1_pattern,R2_pattern=readTrack_naming_convention.split('|')


if 'fastp_param' in config.key() :
    fastp_param=config['fastp_param']
else:
    fastp_param=''

## retrieve name tag list and input path

if os.path.exists(input_): # check for input_ exists, then it’ll check if 
    # input_ is a file or a directory. If it’s a directory, it’ll take all file
    # that end with the file format provided. Otherwise, it’ll take only the file 
    # input_ as input.
    if os.path.isfile(path):  
        print("It is a normal file.")
        print('input_ detected as a single file')
        sample_list=input_
        print(f'input file considered: {input_}')
        sample_name=re.sub('.*/','',sample_list)
        name_tag=sample_name.replace(format_file,'')
        input_path=input_.replace(sample_list,'')
    elif os.path.isdir(path): 
        print("It is a directory.")
        print('input_ detected as a directory, file format provided will be considered')
        sample_list=[os.path.join(input_,f) for f in os.listdir(input_) if f.endwith(format_file) ]
        sample_list_join='\n'.join(sample_list)
        print(f'input files considered: \n{sample_list_joined}')
        sample_name=[f for f in os.listdir(input_) if f.endwith(format_file)]
        name_tag=[f.replace(format_file,'') for f in sample_name]
    else:
        print("Something went wrong !! The provided file is a special file (socket, FIFO, device file), which is not covered by SeqRepeatPuller")
        exit 0
else:
    print("The input path provided do not exist")
    exit 0
## setup workdir
workdir: input_path

## Determine from format file the output list
output_files=['Xtracted_seq.txt','sorted.bam','sorted.bam.bai']
if format_file!='bam' :
    if aligner=='bismark':
        output_files.append('bismark.all.bam')
    elif aligner=='star':
        output_files.append('star.all.bam')

    if format_file!='merged.fastq.gz':
        output_files.append('merged.fastq.gz')
        name_tag=list(set([re.sub(read_naming_convention,'',x) for x in name_tag]))
            
    readTrack_names=[R1_pattern,R2_pattern]
else:
    readTrack_names=''

if format_file is not in ['fq','fastq','bam'] :
    print("format file provided non allowed. Only bam  or fastq and fq are allowed")
    exit 0

onsucces:
    print("SeqRepeatPuller have finished, no error")

rule all:
    inputs=expand('{input_path}/{name_tag}{readTrack_name}.{format_file}',input_path=input_path,readTrack_name=readTrack_names,name_tag=name_tag,format_file=format_file),
    outputs=expand('{work_dir}/{name_tag}_{output_files}',work_dir=input_path,name_tag=name_tag,output_files=output_files)


rule merge:
    input:
        R1='{name_tag}{R1_pattern}.{format_file}',
        R2='{name_tag}{R2_pattern}.{format_file}'
    params:
        fastp_param
    output:
        '{name_tag}_merged.fastq.gz'
    shell:
        """
        a="{format_file}"
        fastp -i {input.R1} -m -I {input.R2} {params} --merged_out {output}

        """

rule mapping_bismark:
    input:
        '{name_tag}_merged.fastq.gz'
    output:
        bam='{name_tag}_bismark.bam',
        tmp_dir=temp(directory('{name_tag}_bismark_tempDir')),
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
        -N 1 \
        -L 15 \
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

rule sort_and_index:
    input:
        '{name_tag}_{aligner}.bam'
    output:
        sortBam=temp('{name_tag}_{aligner}_sorted.bam'),
        index=temp('{name_tag}_{aligner}_sorted.bam.bai')
    shell:
        """
        samtools sort -o {output.sortBam} {input}
        samtools index {output.sortBam}
        """


rule Xtraction:
    input:
        bam='{name_tag}_{aligner}.bam',
        pos_file=pos_file
    output:
        '{name_tag}_Xtracted_seq.txt'
    shell:
        """
        python3 {binpath}/core/Xtractor.py -b {input.bam} -p {input.pos_file} -o {output}  
        """

rule generate_metadata:
    input:
        '{name_tag}_Xtracted_seq.txt'
    output:
        temp('{name_tag}_metadata.json')
    run:
        import json



        
rule add_metadata_and_zip:
    input:
        xtract='{name_tag}_Xtracted_seq.txt',
        metadata='{name_tag}_metadata.json'
    output:
        'processedData/{name_tag}_Xtracted_seq.zip'
    shell:
        """
        zip {output} {input.xtract} {input.metadata}
        """








