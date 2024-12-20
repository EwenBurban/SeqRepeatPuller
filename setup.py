from setuptools import setup, find_packages

setup(
    name="snakemake_wrapper",
    version="1.0",
    packages=find_packages(),
    package_data={
        "snakemake_wrapper": [
            "pipeline_*/*",   # Include pipeline files
            "binpath_*/*",   # Include binpath scripts
        ]
    },
    entry_points={
        "console_scripts": [
            "SeqRepeatPuller=snakemake_wrapper.wrapper_SeqRepeatPuller:main",
            "get_WG_pos=snakemake_wrapper.wrapper_get_WG_pos:main",
        ],
    },
    include_package_data=True,
    install_requires=[],
)

