import os
import subprocess
import sys
import argparse
import datetime


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Wrapper for indexing reference genmoe depending of the aligner.")
    parser.add_argument("-i","--RefSeqFasta",required=True,help="the reference fasta on which samples must be align")
    parser.add_argument("-a","--aligner",required=True,help="the name of the aligner program, could be bowtie2 or bismark")
    parser.add_argument("-o","--output",required=True,help="output file or directory depending of the aligner")
    return parser.parse_args()




def run_process(args):
    # build and execute the command
    print("warning, for bowtie, output expect a file where for bismark output is a directory")
    
    if args.aligner == "bowtie2":
        cmd = ["bowtie2-build", args.RefSeqFasta, args.output]
    elif args.aligner == "bismark":
        cmd = ["bash", "-c", f"mkdir -p {args.output} && cp {args.RefSeqFasta} {args.output}/ && bismark_genome_preparation --bowtie2 {args.output}"]
    else:
        print("Aligner not supported. Only bowtie2 and bismark are supported.")
        sys.exit(1)
    result = subprocess.run(cmd)

    
    if result.returncode == 0:
        print("indexing completed successfully.")
    else:
        print(f"indexing failed")
        sys.exit(result.returncode)

def main():
    arguments=parse_arguments()

    run_process(arguments)


if __name__ == "__main__":
    main()

