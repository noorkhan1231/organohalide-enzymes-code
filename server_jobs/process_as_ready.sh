#!/bin/bash
# Runs fastp -> MicrobeCensus -> Salmon (all 10 families) on each run AS SOON AS its download is finished (links/ncbi_done/<run>),
# while ncbi_download.sh keeps downloading the others. Same tools and settings as overnight_ambient2.sh. Safe to restart.
# Usage (in Noor/Ambient):  nohup bash process_as_ready.sh chicken_topup_core.tsv > ready.log 2>&1 &
cd "$(dirname "$0")"
L=${1:-chicken_topup_core.tsv}; JOBS=${JOBS:-6}; TH=${TH:-5}      # 6 x 5 = 30 threads
E=/home/students/yangkaiwen/.conda/envs
FASTP=$E/transcriptome/bin/fastp; SALMON=$E/salmon_env/bin/salmon; MC=$E/microbecensus_env/bin/run_microbe_census.py
FA=$(cat fastp_args.txt 2>/dev/null || echo "-q 20 -u 20 -l 50 --trim_front1 10 --trim_tail1 0"); SA=$(cat salmon_args.txt)
mkdir -p clean qc mc ambient_quant links/proc_done
tail -n +2 "$L" | cut -f3 | sort -u > links/chicken_new_runs.txt
echo "start $(date); runs: $(wc -l < links/chicken_new_runs.txt); fastp args: $FA"
proc(){ r=$1
  if [ ! -s qc/$r.json ]; then
    rm -f clean/${r}_1.fq.gz clean/${r}_2.fq.gz
    for x in "" ".gz"; do [ -s raw_reads/${r}_2.fastq$x ] && { P1=raw_reads/${r}_1.fastq$x; P2=raw_reads/${r}_2.fastq$x; break; }; done
    for x in "" ".gz"; do [ -s raw_reads/${r}.fastq$x ] && S1=raw_reads/${r}.fastq$x; done
    if [ -n "${P2:-}" ]; then $FASTP -i $P1 -I $P2 -o clean/${r}_1.fq.gz -O clean/${r}_2.fq.gz $FA -w $TH -j qc/$r.json -h qc/$r.html >/dev/null 2>qc/$r.err
    else $FASTP -i $S1 -o clean/${r}_1.fq.gz $FA -w $TH -j qc/$r.json -h qc/$r.html >/dev/null 2>qc/$r.err; fi
    [ -s qc/$r.json ] && rm -f raw_reads/${r}_1.fastq raw_reads/${r}_2.fastq raw_reads/${r}.fastq   # plain fastq no longer needed
  fi
  [ -s qc/$r.json ] || { echo "FAILED fastp $r: $(tail -c 300 qc/$r.err | tr '\n' ' ')"; return; }
  mcpid=""; if [ ! -s mc/$r.txt ]; then $MC -n 2000000 -t 1 clean/${r}_1.fq.gz mc/$r.txt >/dev/null 2>mc/$r.err & mcpid=$!; fi   # MicrobeCensus runs alongside Salmon
  while IFS=$'\t' read fam idx; do
    o=ambient_quant/$fam/$r; [ -s $o/quant.sf ] && continue; mkdir -p ambient_quant/$fam
    if [ -f clean/${r}_2.fq.gz ]; then $SALMON quant -i $idx $SA -1 clean/${r}_1.fq.gz -2 clean/${r}_2.fq.gz -p $TH -o $o >/dev/null 2>$o.err
    else $SALMON quant -i $idx $SA -r clean/${r}_1.fq.gz -p $TH -o $o >/dev/null 2>$o.err; fi
    [ -s $o/quant.sf ] || echo "FAILED salmon $fam $r"
  done < salmon_indices.tsv
  [ -n "$mcpid" ] && wait $mcpid
  [ -s mc/$r.txt ] || echo "FAILED microbecensus $r"
  touch links/proc_done/$r; echo "done $r $(date)"; }
export -f proc; export FASTP SALMON MC FA SA TH
# continuous scheduler: keeps JOBS runs going at all times; a new run starts the moment a slot is free and its download is finished
[ "${KEEP:-0}" = 1 ] || rm -rf links/proc_started; mkdir -p links/proc_started    # KEEP=1: do not forget runs that another (older) scheduler is still processing
while true; do
  started=0
  while read r; do
    [ -e links/ncbi_done/$r ] && [ ! -e links/proc_done/$r ] && [ ! -e links/proc_started/$r ] || continue
    while [ "$(jobs -rp | wc -l)" -ge "$JOBS" ]; do sleep 10; done
    touch links/proc_started/$r; started=1
    bash -c "proc $r" &
  done < links/chicken_new_runs.txt
  if [ $started = 0 ]; then
    grep -q "NCBI runs finished" ncbi_chicken.log 2>/dev/null && [ -z "$(jobs -rp)" ] && break     # downloads over, nothing running
    sleep 30
  fi
done
wait
echo "ALL DONE $(date): processed $(ls links/proc_done | wc -l) of $(wc -l < links/chicken_new_runs.txt); quant dirs $(ls ambient_quant/*/*/quant.sf | wc -l)"
touch CHICKEN_PROC_DONE
