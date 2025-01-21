import os
import re
import sys
#### mandatory arguments
input_=config['input_']
format_file= config['format_file']
pos_file=config['pos_file']
binpath=config['binpath']
work_dir=config['workdir']
git_commit=config['git_commit']
#### optional arguments
if 'genome_folder' in config.keys():
    genome_folder=config['genome_folder']

if 'read_naming' in config.keys() :
    readTrack_naming_convention=config['read_naming']
else : 
    readTrack_naming_convention='_R1|_R2'
R1_pattern,R2_pattern=readTrack_naming_convention.split('|')

if 'aligner' in config.keys():
    aligner=config['aligner']
else :
    aligner=''

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
        # Extract the file name without the path

        # Extract the base name without format (e.g., `.fastq.gz`)
        sample_name = os.path.basename(input_)  # Get the file name with extension
        name_tags = [sample_name.replace('.' + format_file, '')]
        input_path = os.path.dirname(input_) + "/"
        if sample_name.endswith("." + format_file):
            name_tags = [sample_name[: -len("." + format_file)]]
        else:
            print("Warning: The format_file does not match the sample_name as expected.")
            name_tags = [sample_name]  # Fallback to the full sample name

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
output_files=['Xtracted_seq.txt','sorted.bam','sorted.bam.bai',"Xtracted_seq_metadata.json"]
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
    else:
        print('aligner provided is not known. Only bismark, bowtie2 are supported for the moment. If you desire to use another aligner, you can provide only the bam file to SeqRepeatPuller')
        sys.exit(0)
    

    if format_file!='merged.fastq.gz' and format_file!='merged.fq.gz':
        output_files.append('merged.fastq.gz')
        name_tags=list(set([re.sub(read_naming_convention,'',x) for x in name_tags]))
            
    readTrack_names=[R1_pattern,R2_pattern]
    dir1 = input_path
    dir2 = work_dir
    dir_final = work_dir
else:
    readTrack_names=''

if format_file not in ['fq.gz', 'fastq.gz', 'bam','merged.fastq.gz','merged.fq.gz']:
    print("format file provided non allowed. Only bam or fastq and fq are allowed")
    sys.exit(0)

if aligner != '':
    aligner="_" + aligner

wildcard_constraints:
    aligner=aligner,
    name_tag='|'.join(name_tags),

rule all:
    input:
        outputs=expand('{dir_final}/{name_tags}_{output_files}',dir_final=dir_final,name_tags=name_tags,output_files=output_files)
    shell:
        """
        """


rule merge:
    input:
        r1=f'{dir1}/{{name_tag}}{R1_pattern}.{format_file}',
        r2=f'{dir1}/{{name_tag}}{R2_pattern}.{format_file}'
    params:
        fastp_param
    threads: 1
    output:
        f'{dir_final}/{{name_tag}}.merged.fastq.gz'
    shell:
        """
        a="{format_file}"
        fastp -i {input.r1} -m -i {input.r2} {params} --merged_out {output}

        """

rule mapping_bismark:
    input:
        f'{dir_final}/{{name_tag}}.merged.fastq.gz'
    threads: 5
    output:
        bam=f'{dir_final}/{{name_tag}}_bismark.bam',
        tmp_dir=temp(directory(f'{dir_final}/.tmp_bismark_{{name_tag}}/')),
        ambig=temp(f'{dir_final}/{{name_tag}}_bismark.ambig.bam')
    shell:
        """
            mkdir -p {output.tmp_dir} ; bismark \
            --fastq \
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
            mv {output.bam} {output.tmp_dir}/temporary.bam
            samtools merge -o {output.bam} {output.tmp_dir}/temporary.bam {output.ambig}
        """

rule mapping_bowtie2:
    input:
        f'{dir_final}/{{name_tag}}.merged.fastq.gz'
    threads: 1
    output:
        f'{dir_final}/{{name_tag}}_bowtie2.bam'
    shell:
        """
        bowtie2 -x {genome_folder} -U {input} | samtools view -b - > {output}
        """


rule sort_and_index:
    input:
        f'{dir2}/{{name_tag}}{aligner}.bam'
    threads: 1
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
    threads: 1
    output:
        f'{dir_final}/{{name_tag}}_Xtracted_seq.txt'
    shell:
        """
        python3 {binpath}/Xtractor.py -b {input.bam} -p {input.pos_file} -o {output} -N {wildcards.name_tag} 
        """

rule generate_metadata:
    input:
        rules.xtraction.output
    threads: 1
    output:
        f'{dir_final}/{{name_tag}}_Xtracted_seq_metadata.json'
    shell:
        """
        python3 {binpath}/generate_metadata.py -i {input} -c "{config}" -o {output}
        """
