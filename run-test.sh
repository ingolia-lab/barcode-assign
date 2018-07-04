#!/bin/bash

rm -r scratch

mkdir -p scratch

./target/debug/bc-seqs --barcodes data/test-bctrim.1.fq --sequences data/test-bctrim.2.fq --outbase scratch/test

./target/debug/bc-align --reference data/CYH2.fa --barcoded-fastq scratch/test_barcoded.fq --outbase scratch/test

