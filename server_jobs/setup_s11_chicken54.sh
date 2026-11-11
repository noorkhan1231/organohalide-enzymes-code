#!/bin/bash
# Build an inputs folder for s11_ambient.py = the original inputs + the 54 further chicken surroundings runs. The original inputs are not changed.
# usage: bash setup_s11_chicken54.sh <FINAL inputs dir> <new inputs dir> <folder with the chicken54_* count files> <topup_candidates.tsv with the 54 rows> [null_summary_chicken54.tsv]
set -eu
F=$1; D=$2; X=$3; TOPUP=$4; NS=${5:-none}
rm -rf $D; mkdir -p $D/original $D/bundle/Ambient
for p in $F/*; do b=$(basename $p); [ "$b" = original ] || [ "$b" = bundle ] || ln -s $p $D/$b; done
for p in $F/bundle/*; do b=$(basename $p); [ "$b" = Ambient ] || ln -s $p $D/bundle/$b; done
for p in $F/bundle/Ambient/*; do b=$(basename $p); [ "$b" = topup_candidates.tsv ] || ln -s $p $D/bundle/Ambient/$b; done
cp $TOPUP $D/bundle/Ambient/topup_candidates.tsv
for p in $F/original/*; do b=$(basename $p); case $b in per_target_counts.tsv.gz|quant_runs.tsv|genome_equivalents.tsv|null_summary.tsv) ;; *) ln -s $p $D/original/$b;; esac; done
{ zcat $F/original/per_target_counts.tsv.gz; zcat $X/chicken54_counts.tsv.gz | tail -n +2; } | gzip > $D/original/per_target_counts.tsv.gz
{ cat $F/original/quant_runs.tsv; tail -n +2 $X/chicken54_quant_runs.tsv; } > $D/original/quant_runs.tsv
{ cat $F/original/genome_equivalents.tsv; tail -n +2 $X/chicken54_ge.tsv; } > $D/original/genome_equivalents.tsv
if [ "$NS" != none ]; then { cat $F/original/null_summary.tsv; tail -n +2 $NS; } > $D/original/null_summary.tsv
else cp $F/original/null_summary.tsv $D/original/null_summary.tsv; fi
echo "inputs ready: $D (null for the 54 runs: $NS)"
