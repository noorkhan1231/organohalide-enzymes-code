import sys, glob, os, re, pandas as pd, numpy as np
G, sd5, lens, out, mode = sys.argv[1:6]
x = pd.ExcelFile(sd5)
fs = pd.concat([x.parse(s) for s in ['fold_search_flagship', 'fold_search_indicator']])
t = pd.read_csv(sys.argv[6], sep='\t')
cur = set(t.hit_id)
sid2hit = dict(zip(fs.safe_id, fs.original_hit_id))
L = pd.read_csv(lens).drop_duplicates('original_hit_id').set_index('original_hit_id').len
os.makedirs(out, exist_ok=True); rows = []
for f in sorted(glob.glob(f'{G}/*/structural_comparison_results/groups/*/tmscore.tsv')):
    grp = f.split('/')[-2]
    d = pd.read_csv(f, sep='\t', header=None, names=['q', 't', 'aln', 'qtm', 'ttm', 'lddt'])
    d = d[d.q != d.t]
    d['sim'] = (d.qtm + d.ttm) / 2
    ids = sorted(set(d.q) | set(d.t))
    if mode == 'curated': ids = [i for i in ids if sid2hit.get(i) in cur]
    M = pd.DataFrame(np.nan, index=ids, columns=ids)
    for r in d.itertuples():
        if r.q in M.index and r.t in M.index:
            M.loc[r.q, r.t] = r.sim if np.isnan(M.loc[r.q, r.t]) else (M.loc[r.q, r.t] + r.sim) / 2
    M = M.combine_first(M.T)
    A = M.values.copy(); np.fill_diagonal(A, 1.0); M = pd.DataFrame(A, index=ids, columns=ids)
    miss = int(M.isna().sum().sum()); M = M.fillna(0.0)     # pairs without a Foldseek alignment: similarity 0
    D = 1 - M
    host = [re.sub(r'_\d+$', '', i).split('_')[-1] for i in ids]
    meta = pd.DataFrame(dict(id=ids, host=host, length=[L.get(sid2hit.get(i), np.nan) for i in ids]))
    D.to_csv(f'{out}/{grp}_dist.csv'); meta.to_csv(f'{out}/{grp}_meta.csv', index=False)
    rows.append(dict(group=grp, n=len(ids), hosts=meta.host.nunique(), missing_pairs=miss, n_noncurated_removed=sum(sid2hit.get(i) not in cur for i in set(d.q)|set(d.t)) if mode=='curated' else 0))
pd.DataFrame(rows).to_csv(f'{out}/groups.csv', index=False); print(pd.DataFrame(rows).to_string())
