#!/bin/bash
# JOB 17 v2 (metadata only, NOTHING is downloaded; about 5-10 min): poultry-farm SURROUNDING metagenomes (litter, bedding, house dust, manure, farm soil)
# One small ENA query per term with retries (v1 failed because one query was too large). Gut/faecal/cloacal samples are removed (those are host).
# Run:    cd /home/students/yangkaiwen/Noor && sed -i 's/\r$//' job17_chicken_surroundings_search.sh && bash job17_chicken_surroundings_search.sh 2>&1 | tee job17.log
# Upload: Noor/Topup17/project_summary.tsv  review_table.tsv
N=/home/students/yangkaiwen/Noor; O=$N/Topup17; rm -rf $O; mkdir -p $O/q; cd $O
API=https://www.ebi.ac.uk/ena/portal/api/search
F="run_accession,study_accession,sample_accession,scientific_name,sample_title,isolation_source,host,environment_material,library_layout,library_strategy,library_source,instrument_platform,base_count,country,collection_date,study_title"
BASE='library_source="METAGENOMIC" AND library_strategy="WGS" AND (instrument_platform="ILLUMINA" OR instrument_platform="DNBSEQ") AND base_count>=2000000000 AND base_count<=50000000000'
i=0
for t in "poultry litter" "chicken litter" "broiler litter" "poultry house" "chicken house" "chicken coop" "poultry farm" "chicken farm" "broiler house" "poultry manure" "chicken manure" "poultry bedding" "hatchery" "poultry dust" "chicken dust" "layer house" "broiler farm" "poultry environment"; do
  for fld in sample_title isolation_source scientific_name study_title environment_material; do
    i=$((i+1)); ok=0
    for k in 1 2 3 4; do
      curl -sS -m 240 -G "$API" --data-urlencode "result=read_run" --data-urlencode "query=$BASE AND $fld=\"*$t*\"" --data-urlencode "fields=$F" --data-urlencode "format=tsv" --data-urlencode "limit=0" -o q/$i.tsv 2>/dev/null && head -1 q/$i.tsv | grep -q '^run_accession' && { ok=1; break; }
      sleep 5
    done
    echo "$t / $fld: $([ $ok = 1 ] && echo $(($(wc -l < q/$i.tsv)-1)) runs || echo FAILED)"
  done
done
python3 - <<'P'
import pandas as pd,glob
d=pd.concat([pd.read_csv(f,sep='\t',dtype=str) for f in glob.glob('q/*.tsv') if open(f).readline().startswith('run_accession')]).drop_duplicates('run_accession')
txt=d[['sample_title','isolation_source','scientific_name','environment_material','study_title']].fillna('').agg(' '.join,axis=1)
d=d[~txt.str.contains('gut|caecum|cecum|cecal|faec|fecal|feces|cloaca|intestin|rumen|ileum',case=False)]
d.to_csv('review_table.tsv',sep='\t',index=False)
s=d.groupby('study_accession').agg(runs=('run_accession','size'),country=('country',lambda v:'; '.join(sorted(set(v.dropna())))),title=('study_title','first'),examples=('sample_title',lambda v:' | '.join(list(v.fillna(''))[:2]))).reset_index()
s.to_csv('project_summary.tsv',sep='\t',index=False); print(len(d),'runs in',len(s),'projects'); print(s[['study_accession','runs','country','title']].to_string(index=False))
P
