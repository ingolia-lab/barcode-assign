#!/bin/bash -x

READ1="/mnt/ingolialab/FastQ/171009_25x125_MS1/NINI015_S1_L001_R1_001.fastq.gz"
READ2=`echo ${READ1} | sed s/_R1_/_R2_/`
DATADIR="/mnt/ingolialab/ingolia/Cyh2"
CONSTANT='GCGATAAAAG$'

BCTRIMBASE="${DATADIR}/171009_bctrim"
BCTRIM_R1="${BCTRIMBASE}.1.fq"
BCTRIM_R2="${BCTRIMBASE}.2.fq"

if [[ ! -e ${BCTRIM_R1} ]];
then
    cutadapt -a ${CONSTANT} \
	   -o ${BCTRIM_R1} -p ${BCTRIM_R2} \
	   --untrimmed-output "${BCTRIMBASE}-reject.1.fq" \
	   --untrimmed-paired-output "${BCTRIMBASE}-reject.2.fq" \
	   ${READ1} ${READ2} > "${BCTRIMBASE}-report.txt"
else
    echo "Skipping barcode trim because ${BCTRIM_R1} exists"
fi
