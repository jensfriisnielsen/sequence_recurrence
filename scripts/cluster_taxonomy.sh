#!/bin/bash

INFILE=$1

OUTDIR=results/taxonomy
mkdir -p $OUTDIR
OUTSTEM=$OUTDIR/taxonomy
time scripts/cluster_taxonomy.py $INFILE data/blastn+x_besthit.txt > ${OUTSTEM}.txt 2> ${OUTSTEM}.time
