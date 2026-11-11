#!/bin/bash
# System-matched null genes (Null_system_matched/null_index, the same index used for the final analysis) for the 54 new chicken surroundings runs.
# Same settings as job2 Part A: first 10 M read pairs, Salmon -l A. Then writes the summary rows in the null_summary.tsv format.
# Run in Noor/Ambient:  sed -i 's/\r$//' null_sysmatched_chicken54.sh && nohup bash null_sysmatched_chicken54.sh > nullsm54.log 2>&1 &
set -u
N=/home/students/yangkaiwen/Noor; SALMON=/home/students/yangkaiwen/.conda/envs/salmon_env/bin/salmon
cd $N/Ambient; W=$N/Null_system_matched
[ -d $W/null_index ] || { echo "MISSING $W/null_index"; exit 1; }
NR=10000000; JOBS=9; TH=4; mkdir -p $W/quant/ambient $W/tmp
q(){ r=$1; o=$W/quant/ambient/$r; [ -s $o/quant.sf ] && return; t=$W/tmp/$r; mkdir -p $t
  f1=$N/Ambient/clean/${r}_1.fq.gz; f2=$N/Ambient/clean/${r}_2.fq.gz
  zcat $f1 | head -n $((4*NR)) > $t/1.fq
  if [ -f "$f2" ]; then zcat $f2 | head -n $((4*NR)) > $t/2.fq; $SALMON quant -i $W/null_index -l A -1 $t/1.fq -2 $t/2.fq -p $TH -o $o >/dev/null 2>$o.err
  else $SALMON quant -i $W/null_index -l A -r $t/1.fq -p $TH -o $o >/dev/null 2>$o.err; fi
  rm -rf $t; [ -s $o/quant.sf ] || echo "FAILED $r"; }
export -f q; export SALMON W N NR TH
xargs -a chicken54_ids.txt -P $JOBS -I{} bash -c 'q {}'
echo "quant done $(date): $(for r in $(cat chicken54_ids.txt); do ls $W/quant/ambient/$r/quant.sf; done 2>/dev/null | wc -l) of $(wc -l < chicken54_ids.txt)"
python3 - <<'PY'
import json, os
W='/home/students/yangkaiwen/Noor/Null_system_matched'
ids=[l.strip() for l in open('chicken54_ids.txt') if l.strip()]
with open('null_summary_chicken54.tsv','w') as o:
    o.write('group\trun\tnum_processed_subsample\tgene_system\treads\trpk\tn_genes\tn_detected\n')
    for r in ids:
        p=f'{W}/quant/ambient/{r}/quant.sf'
        if not os.path.exists(p): continue
        try: npro=json.load(open(os.path.dirname(p)+'/aux_info/meta_info.json'))['num_processed']
        except Exception: npro=''
        s={}
        with open(p) as fh:
            next(fh)
            for l in fh:
                c=l.split('\t'); k=c[0].split('|')[0]; L=float(c[1]); n=float(c[4])
                a=s.setdefault(k,[0,0.0,0,0]); a[0]+=n; a[1]+=n/(L/1000); a[2]+=1; a[3]+=n>0
        for k,v in s.items(): o.write('\t'.join(map(str,('ambient',r,npro,k,*v)))+'\n')
print('written null_summary_chicken54.tsv')
PY
wc -l null_summary_chicken54.tsv; head -3 null_summary_chicken54.tsv
