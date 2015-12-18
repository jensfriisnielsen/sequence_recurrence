#!/bin/bash

INCLSTR=$1
BNAME=$(basename $1)
OUTDIR=results/feature_associations
mkdir -p $OUTDIR
INBAK=${INCLSTR%.clstr}.bak.clstr
FEATURES=data/features

# collapse clusters into table
OUTCLSTR=$OUTDIR/$(basename $INCLSTR)
DATASETS=$OUTCLSTR.datasets
cat $INCLSTR | sed -r 's/ +/\t/g' | cut -f1,3 | tr -d '.>' | cut -f1 -d_ | sort -k1,1g | scripts/join_lines.py > $DATASETS

# create fisher exact test association matrix
# cluster0 f001 c+f+, c-f+, c+f-, c-f-
MAT=${DATASETS}.feature_association_matrix
cat $DATASETS | scripts/cluster_feature_associations.py --all-datasets data/all_datasets $FEATURES/f???.{+,-} > $MAT

# find pvalues from feature association matrix
PVA=$MAT.pvalue
R --slave --vanilla --args $MAT $PVA < scripts/cluster_feature_associations.R
mv $PVA $OUTDIR/pval.txt
