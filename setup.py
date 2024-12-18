from setuptools import setup, find_packages

setup(
    name="snakemake_wrapper",
    version="1.0",
    packages=find_packages(),
    entry_points={
        "console_scripts": [
            "wrapper_tool1=snakemake_wrapper.wrapper_SeqRepeatPuller:main",
            "wrapper_tool2=snakemake_wrapper.wrapper_get_WG_pos:main",
        ],
    },
    include_package_data=True,
    install_requires=[],
)

