#!/bin/bash
# Two-stage NCBI download: stage 1 = prefetch, PJOBS at a time (network only, tiny disk use);
#                          stage 2 = fasterq-dump, only CONV at a time, temp files in RAM (/dev/shm) when >= 40 GB is free there, else on disk.
# Output: plain fastq in raw_reads (no gzip; process_as_ready.sh reads it and deletes it after fastp). Safe to restart.
# Usage (in Noor/Ambient):  PJOBS=30 CONV=4 nohup bash ncbi_two_stage.sh chicken_topup_core.tsv > ncbi_chicken.log 2>&1 &
cd "$(dirname "$0")"; mkdir -p raw_reads links/ncbi_done links/sra_done sra ncbi_tmp
TK=/home/students/yangkaiwen/sratoolkit.2.10.8-ubuntu64/bin
PJ=${PJOBS:-30}; CV=${CONV:-4}; ET=${ETHREADS:-4}; RAMNEED=${RAMGB:-40}
L=${1:-chicken_topup_core.tsv}
col=$(head -1 "$L" | tr '\t' '\n' | grep -nx -E 'Run|run_accession' | head -1 | cut -d: -f1)
tail -n +2 "$L" | cut -f"$col" | sort -u > links/ncbi_runs.txt; echo "NCBI runs: $(wc -l < links/ncbi_runs.txt)"
dl(){ r=$1
  [ -e links/ncbi_done/$r ] && return; [ -e links/sra_done/$r ] && return
  for k in 1 2 3 4 5; do $TK/prefetch $r -O sra --max-size u >/dev/null 2>&1 && [ -s sra/$r/$r.sra ] && break; sleep 30; done
  [ -s sra/$r/$r.sra ] && touch links/sra_done/$r || echo "FAILED prefetch $r"; }
export -f dl; export TK
xargs -a links/ncbi_runs.txt -P $PJ -I{} bash -c 'dl {}' & DLPID=$!
conv(){ r=$1
  freegb=$(df --output=avail -k /dev/shm | tail -1 | awk '{print int($1/1048576)}')
  if [ "$freegb" -ge "$RAMNEED" ]; then t=/dev/shm/$USER/$r; else t=$PWD/ncbi_tmp/$r; fi
  mkdir -p $t; rm -f raw_reads/${r}*.fastq*
  $TK/fasterq-dump sra/$r/$r.sra -e $ET --split-files -t $t -O raw_reads >/dev/null 2>&1 \
    || { echo "FAILED convert $r"; rm -rf $t raw_reads/${r}*.fastq*; touch links/conv_failed_$r; return; }
  rm -rf $t sra/$r; touch links/ncbi_done/$r; echo "converted $r ($([ "$t" = "/dev/shm/$USER/$r" ] && echo RAM || echo disk) temp) $(date +%H:%M)"; }
export -f conv; export TK ET RAMNEED
mkdir -p links/conv_started; rm -f links/conv_started/* links/conv_failed_*
while true; do
  started=0
  while read r; do
    [ -e links/sra_done/$r ] && [ ! -e links/ncbi_done/$r ] && [ ! -e links/conv_started/$r ] || continue
    while [ "$(jobs -rp | wc -l)" -ge $((CV+1)) ]; do sleep 10; done          # +1 = the download pipeline itself
    touch links/conv_started/$r; started=1; bash -c "conv $r" &
  done < links/ncbi_runs.txt
  if [ $started = 0 ]; then
    kill -0 $DLPID 2>/dev/null || { [ "$(jobs -rp | wc -l)" -le 0 ] && break; }
    sleep 30
  fi
done
wait
echo "NCBI runs finished: $(ls links/ncbi_done | wc -l) of $(wc -l < links/ncbi_runs.txt) (list includes earlier runs)"
touch links/NCBI_DONE; echo done
