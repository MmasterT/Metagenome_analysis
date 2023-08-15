""" Installation of """
#!/usr/bin/env python

import glob
from distutils.util import convert_path
from setuptools import setup




requirements = [line.rstrip() for line in open("requirements.txt", "rt")]

main_ns = {}
ver_path = convert_path('yaragenome/version.py')
with open(ver_path, encoding="utf-8") as ver_file:
    exec(ver_file.read(), main_ns)

setup(
    name="eihic",
    version=main_ns['__version__'],
    description=" Metagenome assembled genomes pipeline",
    author="Mariano Olivera Fedi",
    author_email="Mariano.olivera-fedi@aerlham.ac.uk, mariano.olivera.fedi@hotmail.com",
    url="https://github.com/MmasterT/Metagenome_analysis",
    license="CC BY-NC 4.0",
    zip_safe=False,
    keywords="snakemake MAGs Illumina ONT PacBio assembler",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Topic :: Scientific Engineering :: Bio/Informatics",
        "License :: OSI Approved :: Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3.9",
    ],
    entry_points={"console_scripts": ["eihic = yaragenome.__main__:main"]},
    install_requires=requirements,
    packages=["yaragenome", "yaragenome.etc", "yaragenome.scripts", "yaragenome.workflows"],
    scripts=[script for script in glob.glob("yaragenome/scripts/*")],
    package_data={
        "yaragenome.workflows": ["assembler_short_reads.smk",
            "assembler_short_reads.smk",
            "binning_short_reads"],
        "yaragenome.etc": ["hpc_config.json", 
            "run_config.yaml", 
            "samples.csv"],
    },
    include_package_data=True,
)
