"""Step 3. Integrity check of every metagenome result set (original + jobs 6, 6b, 14, 16) and the metatranscriptomes.
usage: python3 s01_integrity.py <inputs_dir> <out_dir>
Checks: 10 families per run; n_targets per family = index size; num_processed identical across the 10 families of a run;
genome equivalents present and > 0; no run in two sets; selection-table host = count host; counts only on known targets.
"""
import sys, os, gzip, pandas as pd, numpy as np
I, O = sys.argv[1:3]; os.makedirs(O, exist_ok=True)
IDX = dict(RDase=276, HyDase=688, GSTDase=3209, OxDase=28, DehyDase=5, MetDase=2149, FlaHase=2052, VanHase=724, NHFeHase=37, SAMHase=19)
rep = []
def say(*a):
    s = ' '.join(str(x) for x in a); print(s); rep.append(s)

sets = {}
q = pd.read_csv(f'{I}/original/quant_runs.tsv', sep='\t'); q = q[q.dataset == 'host']
g = pd.read_csv(f'{I}/original/genome_equivalents.tsv', sep='\t'); g = g[g.dataset == 'host']
sets['original'] = (q, g, f'{I}/original/per_target_counts.tsv.gz')
for j in ['job6', 'job6b', 'job14', 'job16']:
    if not os.path.exists(f'{I}/{j}/extra_quant_runs.tsv'):
        say(f'[{j}] MISSING result files'); continue
    sets[j] = (pd.read_csv(f'{I}/{j}/extra_quant_runs.tsv', sep='\t'), pd.read_csv(f'{I}/{j}/extra_ge.tsv', sep='\t'), f'{I}/{j}/extra_counts.tsv.gz')

meta = pd.read_csv(f'{I}/original/target_meta.tsv', sep='\t')
known = set(zip(meta.family, meta.target))
allruns = []
for name, (q, g, cf) in sets.items():
    runs = q.drop_duplicates('run')[['run', 'host_dir']]
    say(f'\n[{name}] runs {len(runs)} by host {runs.host_dir.value_counts().sort_index().to_dict()}')
    nf = q.groupby('run').family.nunique(); say('  runs without all 10 families:', int((nf != 10).sum()))
    bad = q[q.n_targets != q.family.map(IDX)]; say('  family rows with wrong n_targets:', len(bad), bad.family.value_counts().to_dict() if len(bad) else '')
    npv = q.groupby('run').num_processed.nunique(); say('  runs whose num_processed differs between families:', int((npv > 1).sum()))
    say('  runs with num_processed missing/0:', int(q.groupby('run').num_processed.max().fillna(0).eq(0).sum()))
    ge = g.drop_duplicates('run').set_index('run').genome_equivalents
    ge = pd.to_numeric(ge, errors='coerce')
    miss = runs[~runs.run.isin(ge[ge > 0].index)]
    say('  runs without genome equivalents > 0:', len(miss), list(miss.run)[:10])
    say('  duplicated run rows (run x family):', int(q.duplicated(['run', 'family']).sum()))
    c = pd.read_csv(cf, sep='\t', usecols=['dataset', 'family', 'run', 'target', 'NumReads'])
    c = c[c.dataset.isin(['host', 'extra'])]
    unk = c[[k not in known for k in zip(c.family, c.target)]]
    say(f'  count rows {len(c):,}; targets unknown to target_meta: {unk.target.nunique()} ({unk.NumReads.sum() / max(c.NumReads.sum(), 1) * 100:.3f}% of reads)')
    say('  count runs not in quant table:', len(set(c.run) - set(runs.run)))
    runs = runs.assign(set=name, num_processed=runs.run.map(q.groupby('run').num_processed.max()), GE=runs.run.map(ge))
    allruns.append(runs)

A = pd.concat(allruns)
d = A[A.run.duplicated(keep=False)].sort_values('run')
say(f'\nruns appearing in more than one set: {d.run.nunique()}')
if len(d): say(d.to_string(index=False))
A.to_csv(f'{O}/all_runs_by_set.tsv', sep='\t', index=False)
say('\nTOTAL unique runs', A.run.nunique(), A.drop_duplicates('run').host_dir.value_counts().sort_index().to_dict())

# selection tables vs counted runs
for j, f in [('job6', 'selected_all.tsv'), ('job14', 'final_selection.tsv'), ('job16', 'final_selection.tsv')]:
    if j not in sets: continue
    s = pd.read_csv(f'{I}/{j}/{f}', sep='\t')
    rc = 'run_accession' if 'run_accession' in s.columns else s.columns[2]
    hc = 'host_system' if 'host_system' in s.columns else s.columns[0]
    r = A[A.set == j].set_index('run').host_dir
    s = s.drop_duplicates(rc)
    hs = s.set_index(rc)[hc].reindex(r.index)
    say(f'[{j}] selection rows {len(s)} | counted {len(r)} | selected but not counted {int((~s[rc].isin(r.index)).sum())} | counted but not selected {int(hs.isna().sum())} | host mismatches {int((hs.notna() & (hs != r)).sum())}')
    nc = s[~s[rc].isin(r.index)]
    if len(nc): say('   not counted:', nc[[hc, rc]].values.tolist()[:20])

# metatranscriptomes
mq = pd.read_csv(f'{I}/metatx/metatx_quant_runs.tsv', sep='\t')
say(f'\n[metatx] runs {mq.run.nunique()} | runs without 10 families {int((mq.groupby("run").family.nunique() != 10).sum())} | wrong n_targets {int((mq.n_targets != mq.family.map(IDX)).sum())} | num_processed differs {int((mq.groupby("run").num_processed.nunique() > 1).sum())}')
say('  host label:', mq.drop_duplicates('run').host_dir.value_counts(dropna=False).to_dict())
open(f'{O}/integrity_report.txt', 'w').write('\n'.join(rep) + '\n')
