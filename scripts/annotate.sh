#!/bin/bash

module load use.own
module load pr_46500
module load ngs anaconda
module load data_trace

INFILE=$1
ASSOC_FILE=$2

< $ASSOC_FILE cat \
    | cut -f1,2,7 \
    | scripts/annotate_cluster_with_taxonomy_from_differential_fractions.py $INFILE 
