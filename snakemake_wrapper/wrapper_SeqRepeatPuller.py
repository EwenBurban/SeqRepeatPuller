import os
import subprocess
import sys
import argparse
import datetime
def get_paths():
    """Resolve paths for SeqRepeatPuller pipeline and binpath."""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    pipeline_dir = os.path.join(script_dir, "pipeline_SeqRepeatPuller")
    binpath_dir = os.path.join(script_dir, "binpath_SeqRepeatPuller")
    snakefile = os.path.join(pipeline_dir, "SeqRepeatPuller.py")

    if not os.path.exists(snakefile):
        print("Error: Pipeline files for SeqRepeatPuller not found.")
        sys.exit(1)

    return snakefile, binpath_dir


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Wrapper for running a Snakemake pipeline.")
    parser.add_argument(
        "-c", "--config", type=str, required=True, help="Path to the Snakemake configuration file (e.g., config.yaml)."
    )
    parser.add_argument( # TODO : to remove as input path is the workdir
        "-w", "--working-dir", type=str, default=os.getcwd(), help="Working directory for the Snakemake pipeline."
    )
    parser.add_argument(
        "-j", "--cores", type=int, default=1, help="Number of CPU cores to use."
    )
    parser.add_argument(
        "-n", "--dry-run", action="store_true", help="Perform a dry-run to test the pipeline without executing."
    )
    parser.add_argument(
        "-k", "--keep-going", action="store_true", help="Continue executing independent jobs even if one fails."
    )
    parser.add_argument(
        "-p", "--profile", type=str, help="Snakemake profile to use (e.g., for cluster execution)."
    )
    parser.add_argument(
        "--unlock", action="store_true", help="Unlock the working directory in case of a previous failure."
    )
    # TODO : add an argparse to put the sbatch command + test it on IFB
    return parser.parse_args()


def run_snakemake(args,binpath,snakefile):
    """Build and execute the Snakemake command."""
    cmd = ["snakemake"]

    # Add Snakefile
    cmd += ["--snakefile", snakefile]

    # Add configuration file
    cmd += ["--configfile", args.config]

    # Add cores
    cmd += ["--cores", str(args.cores)]

    # Add working directory
    if args.working_dir:
        cmd += ["--directory",str(args.working_dir)]
    # Add profile if specified
    if args.profile:
        cmd += ["--profile", args.profile]

    # Handle special flags
    if args.dry_run:
        cmd += ["--dry-run"]
    if args.keep_going:
        cmd += ["--keep-going"]
    if args.unlock:
        cmd += ["--unlock"]

    cmd += ["--config", f"binpath={binpath}"]
    cmd += ["--config", f"workdir={args.working_dir}"]

    # Add logs for better debugging
    # Get current time details
    now = datetime.datetime.now()
    month=now.month
    year=now.year
    day = now.day
    hour = now.hour
    minute = now.minute
    second = now.second
    log_file = os.path.join(args.working_dir, f"SeqRepeatPuller_{year:02d}{month:02d}{day:02d}{hour:02d}{minute:02d}{second:02d}.log")
    print(f"Logging output to {log_file}")

    # Run Snakemake
    with open(log_file, "w") as log:
        result = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
    
    if result.returncode == 0:
        print("Snakemake pipeline completed successfully.")
    else:
        print(f"Snakemake pipeline failed. Check the log file at {log_file} for details.")
        sys.exit(result.returncode)

def main():
    snakefile, binpath = get_paths()

    arguments=parse_arguments()

    run_snakemake(arguments,binpath,snakefile)

    print("Running Tool 1 pipeline...")

if __name__ == "__main__":
    main()

