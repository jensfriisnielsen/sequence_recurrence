#!/bin/bash

# -c alignment percent-id
# -G 1 global alignment
# -M 0 unlimited memory usage (assigned by queueing system)
# -T threads
# -aS minimum alignment-coverage of shorter sequence
# -B 0 sequences are stored in RAM
# -p 1 print overlap in .clstr
# -g 1 slow accurate mode (most similar cluster)
# -r 1 both +/+ and +/- strand alignments
OUTDIR=results/clustering
mkdir -p $OUTDIR
PARF=$(mktemp)
cat << EOF | tr -d '\t -.' > $PARF
$@
EOF
PAR=$(cat $PARF)
rm -f $PARF
OFILE=$OUTDIR/contigs_all.cdhit.$PAR
if [ -f $OFILE.exit ]; then
    EXIT=$(cat $OFILE.exit)
    exit $EXIT
fi
cd-hit-est -i data/contigs_all.fsa -o $OFILE $@ 2> $OFILE.stderr  > $OFILE.stdout
echo $? > $OFILE.exit
