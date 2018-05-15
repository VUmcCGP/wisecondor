#!/usr/bin/env bash

# PARAMETERS
CORES=6
INPUT_DIR="path/to/convert.npz" # existing (non-empty) folder, containing reference .npz files
OUTPUT_DIR="path/to/newref.npz" # existing (empty) folder
REF_SIZES="15 50 100 200 500 1000" # space separated list (kb)
RELEASE="hg38" # reference used to create bam files (solely used for reference filename)
GENDER="F" # gender of the of the cases at NPZ_INPUT_DIR

# SCRIPT

for REF in ${REF_SIZES}
do
    echo "Creating reference at bins size ${REF} kb"

    WisecondorX newref ${INPUT_DIR}/*.npz \
    ${OUTPUT_DIR}/reference.${RELEASE}.${GENDER}.${REF}kb.npz \
    -binsize ${REF}000 -cpus ${CORES} -gender ${GENDER}
done