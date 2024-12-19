import os
import re
import sys
#### mandatory arguments
input_=config['input_']
format_file= config['format_file']
pos_file=config['pos_file']
binpath=config['binpath']
work_dir=config['workdir']

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
        name_tags = [sample_name.replace('.' + format_file,'')]
        input_path = input_.replace(name_tags[0]+'.' + format_file,'')
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
output_files=['Xtracted_seq.txt','sorted.bam','sorted.bam.bai',"Xtracted_seq.zip"]
dir1=input_path
dir2=input_path
dir_final=work_dir
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
    dir2 = work_dir
    dir_final = work_dir
else:
    readTrack_names=''

if format_file not in ['fq', 'fastq', 'bam']:
    print("format file provided non allowed. Only bam or fastq and fq are allowed")
    sys.exit(0)


wildcard_constraints:
    aligner=aligner,
    name_tag='|'.join(name_tags),

rule all:
    input:
 #       inputs=expand('{input_path}/{name_tag}{readtrack_name}.{format_file}',input_path=input_path,readtrack_name=readtrack_names,name_tag=name_tag,format_file=format_file),
        outputs=expand('{dir_final}/{name_tags}_{output_files}',dir_final=dir_final,name_tags=name_tags,output_files=output_files)
    shell:
        """
        """


rule merge:
    input:
        r1=f'{dir1}/{{name_tag}}{{r1_pattern}}.{{format_file}}',
        r2=f'{dir1}/{{name_tag}}{{r2_pattern}}.{{format_file}}'
    params:
        fastp_param
    output:
        f'{dir_final}/{{name_tag}}_merged.fastq.gz'
    shell:
        """
        a="{format_file}"
        fastp -i {input.r1} -m -i {input.r2} {params} --merged_out {output}

        """
rule mapping:
    input:
        rules.merge.output
    output:
        '{dir_final}/{name_tag}.bam'
    shell:
        """
        if [[ "{aligner}" == "bowtie2" ]]; then
            bowtie2 -x {genome_folder} -U {input} | samtools view -bs4 - > {output}

        elif [[ "{aligner}" == "bismark"; then
            shell("mkdir -p {dir_final}/.tmp_{wildcards.name_tag}/ ; bismark \
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

        elif [[ "{aligner}" == "star"]]; then
            touch {output}'
        """

rule sort_and_index:
    input:
        f'{dir2}/{{name_tag}}.bam'
    output:
        sortbam=temp(f'{dir_final}/{{name_tag}}_sorted.bam'),
        index=temp(f'{dir_final}/{{name_tag}}_sorted.bam.bai')
    shell:
        """
        samtools sort -o {output.sortbam} {input}
        samtools index {output.sortbam}
        """


rule xtraction:
    input:
        bam=rules.sort_and_index.output.sortbam,
        index=rules.sort_and_index.output.index,
        pos_file=pos_file
    output:
        temp(f'{dir_final}/{{name_tag}}_Xtracted_seq.txt')
    shell:
        """
        python3 {binpath}/Xtractor.py -b {input.bam} -p {input.pos_file} -o {output} -N {wildcards.name_tag} 
        """

rule generate_metadata:
    input:
        rules.xtraction.output
    output:
        temp(f'{dir_final}/{{name_tag}}_metadata.json')
    run:
        import json
        from datetime import datetime
        if format_file == "bam":
            aligner="None"

        metadata ={
                "file_name" : input[0],
                "created_at": datetime.now().isoformat(),
                "aligner" : aligner,
                "pos_file": pos_file}
        with open(output[0],"w") as out_file:
            json.dump(metadata,out_file,indent=4)


        
rule add_metadata_and_zip:
    input:
        xtract=rules.xtraction.output,
        metadata=rules.generate_metadata.output
    output:
        f'{dir_final}/{{name_tag}}_Xtracted_seq.zip'
    shell:
        """
        zip {output} {input.xtract} {input.metadata}
        """








