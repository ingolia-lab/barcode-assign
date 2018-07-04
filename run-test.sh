#!/bin/bash

READ1="data/test_R1_001.fastq"
READ2="data/test_R2_001.fastq"
REFERENCE="data/CYH2.fa"
DATADIR="./scratch/"
NAME="test"

READ1_CONSTANT='GCGATAAAAG$'
# READ2_CONSTANT='CCCAACGCTTTTATCGC'

rm -r ${DATADIR}

mkdir -p ${DATADIR}

BCTRIMBASE="${DATADIR}/${NAME}-bctrim"
BCTRIM_R1="${BCTRIMBASE}.1.fq"
BCTRIM_R2="${BCTRIMBASE}.2.fq"

cutadapt -a "${READ1_CONSTANT}" \
         -o "${BCTRIM_R1}" -p "${BCTRIM_R2}" \
         --untrimmed-output "${BCTRIMBASE}-reject.1.fq" \
         --untrimmed-paired-output "${BCTRIMBASE}-reject.2.fq" \
         "${READ1}" "${READ2}" > "${BCTRIMBASE}-cutadapt.txt"

./target/debug/bc-seqs --barcodes "${BCTRIM_R1}" --sequences "${BCTRIM_R2}" --outbase "${DATADIR}/${NAME}"

./target/debug/bc-align --reference "${REFERENCE}" --barcoded-fastq "${DATADIR}/${NAME}_barcoded.fq" --outbase "${DATADIR}/${NAME}"

