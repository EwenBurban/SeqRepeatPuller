from setuptools import setup, find_packages
import subprocess
# Récupérer le nom du commit Git
def get_git_commit():
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], text=True
        ).strip()
    except Exception:
        return "unknown"

# Write the Git commit hash to a file in the package
commit=get_git_commit()
with open("snakemake_wrapper/version.py", "w") as f:
    f.write(f'__commit__= "{commit}"\n')

setup(
    name="snakemake_wrapper",
    version="0.2",
    packages=find_packages(),
    package_data={
        "snakemake_wrapper": [
            "pipeline_*/*",   # Include pipeline files
            "binpath_*/*",   # Include binpath scripts
            "git_commit.txt",# include file containing git version

        ]
    },
    entry_points={
        "console_scripts": [
            "SeqRepeatPuller=snakemake_wrapper.wrapper_SeqRepeatPuller:main",
            "get_WG_pos=snakemake_wrapper.wrapper_get_WG_pos:main",
            "index=snakemake_wrapper.wrapper_index_ref_genome:main",
        ],
    },
    include_package_data=True,
    install_requires=[],
)

