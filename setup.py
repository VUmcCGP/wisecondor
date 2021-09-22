#!/usr/bin/env python

from setuptools import find_packages, setup

version = "1.2.5-dev"
dl_version = "master" if "dev" in version else "{}".format(version)

setup(
    name="WisecondorX",
    version=version,
    author="Matthias De Smet, Lennart Raman",
    author_email="Lennart.raman@ugent.be",
    description="WisecondorX -- an evolved WISECONDOR",
    long_description=__doc__,
    keywords=[
        "bioinformatics",
        "biology",
        "sequencing",
        "NGS",
        "next generation sequencing",
        "CNV",
        "SWGS",
        "Shallow Whole Genome Sequencing",
    ],
    download_url="https://github.com/CenterForMedicalGeneticsGhent/WisecondorX/archive/v{}.tar.gz".format(
        dl_version
    ),
    license="Attribution-NonCommercial-ShareAlike CC BY-NC-SA",
    packages=find_packages("."),
    python_requires=">=3.6",
    include_package_data=True,
    zip_safe=False,
    install_requires=["scipy", "scikit-learn", "pysam", "numpy", "typer"],
    entry_points={"console_scripts": ["WisecondorX = wisecondorX.main:main"]},
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX",
        "Operating System :: Unix",
        "Programming Language :: Python :: 3 :: only",
        "Programming Language :: Python :: 3.6",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
