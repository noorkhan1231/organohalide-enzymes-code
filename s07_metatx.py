"""Step 10. Metatranscriptome expression of the ten families: RPKM = reads per kilobase of target per million reads processed.
usage: python3 s07_metatx.py <inputs_dir> <out_dir>
Run set: the 191 earlier runs (metatx_selected_original.tsv) + 43 added runs (approved_metatx_new_v3.tsv), all processed
identically (first 30 million reads, fastp, Salmon). Library rule (RNA): library source TRANSCRIPTOMIC or METATRANSCRIPTOMIC
with strategy RNA-Seq or OTHER; AMPLICON, miRNA-Seq, WGA and WGS libraries are excluded. Host labels come from the
selection tables (the host_dir column of metatx_quant_runs.tsv is not used, E20). No MicrobeCensus for RNA.
Tests: Kruskal-Wallis per family at run level and study level (study = BioProject or merged series), detection chi-square.
"""
import sys, os, json, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stats_utils import *
from scipy.stats import chi2_contingency
I, O = sys.argv[1:3]; os.makedirs(O, exist_ok=True)
old = pd.read_csv(f'{I}/metatx/metatx_selected_original.tsv', sep='\t')[['host_system', 'study_accession', 'run_accession']]
srch = pd.read_csv(f'{I}/metatx/metatx_all_runs_search.tsv', sep='\t', low_memory=False).drop_duplicates('run_accession')
old = old.merge(srch[['run_accession', 'library_strategy', 'library_source']], on='run_accession', how='left')
new = pd.read_csv(f'{I}/metatx/approved_metatx_new_v3.tsv', sep='\t')[['host_system', 'study_accession', 'run_accession', 'library_strategy']]
new['library_source'] = np.nan
S = pd.concat([old.assign(origin='earlier'), new.assign(origin='added')], ignore_index=True).drop_duplicates('run_accession')
ena_f = f'{I}/metadata/ena_metatx.tsv'
if os.path.exists(ena_f):                                   # ENA record wins where present
    e = pd.read_csv(ena_f, sep='\t', dtype=str).drop_duplicates('run_accession').set_index('run_accession')
    for col in ['library_strategy', 'library_source']:
        S[col] = S.run_accession.map(e[col]).fillna(S[col])
    S['study_accession'] = S.run_accession.map(e.study_accession).fillna(S.study_accession)
ok_src = S.library_source.isna() | S.library_source.isin(['TRANSCRIPTOMIC', 'METATRANSCRIPTOMIC', 'TRANSCRIPTOMIC SINGLE CELL'])
S['include'] = S.library_strategy.isin(['RNA-Seq', 'OTHER']) & ok_src
S['exclude_reason'] = np.where(S.include, '', 'library ' + S.library_strategy.fillna('?') + ' / ' + S.library_source.fillna('?'))
SG = {'PRJNA1016881': 'SG_UCDavis_rhizo_metatx'}
for k in range(1016881, 1016893): SG[f'PRJNA{k}'] = 'SG_UCDavis_rhizo_metatx'
S['study'] = S.study_accession.map(SG).fillna(S.study_accession)
# sample-type decisions (same rules as the metagenomes; metatx_qc_decisions.tsv); 'pending' studies stay out until checked
DQ = pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'metatx_qc_decisions.tsv'), sep='\t')
for _, d in DQ[DQ.decision.isin(['exclude', 'pending'])].iterrows():
    m = (S.run_accession == d.id) if d.level == 'run' else (S.study_accession == d.id)
    if not m.any(): raise SystemExit('metatx decision matched no run: ' + d.id)
    m = m & S.include
    S.loc[m, 'include'] = False; S.loc[m, 'exclude_reason'] = d.code + ': ' + d.reason
S.to_csv(f'{O}/metatx_runs_qc.tsv', sep='\t', index=False)
S = S[S.include]
c = pd.read_csv(f'{I}/metatx/metatx_counts.tsv.gz', sep='\t'); c = c[c.dataset == 'metatx']
q = pd.read_csv(f'{I}/metatx/metatx_quant_runs.tsv', sep='\t'); q = q[q.dataset == 'metatx'].groupby('run').num_processed.max()
c = c[c.run.isin(S.run_accession)].copy()
c['rpk'] = c.NumReads / (c.Length / 1000)
w = c.groupby(['run', 'family']).rpk.sum().unstack(fill_value=0).reindex(columns=FAM, fill_value=0)
w = S.set_index('run_accession')[['host_system', 'study', 'origin']].join(w).fillna({f: 0 for f in FAM})
w['million'] = w.index.map(q) / 1e6
miss = w[w.million.isna() | (w.million <= 0)]
w = w[w.million.notna() & (w.million > 0)]
for f in FAM: w[f] = w[f] / w.million
R = dict(n_runs=int(len(w)), runs_by_host=w.host_system.value_counts().reindex(HOSTS).to_dict(), runs_without_quant=miss.index.tolist(),
         studies_by_host=w.groupby('host_system').study.nunique().reindex(HOSTS).to_dict(), excluded_library=int((~pd.concat([S]).include).sum()))
R['excluded_library'] = int(len(pd.read_csv(f'{O}/metatx_runs_qc.tsv', sep='\t').query('include == False')))
st = w.groupby(['host_system', 'study'])[FAM].mean().reset_index()
g = w.host_system.values; rows = []
for f in FAM:
    H, p = kw(w[f].values, g); Hs, ps = kw(st[f].values, st.host_system.values)
    det = pd.crosstab(w.host_system, w[f] > 0).reindex(HOSTS).fillna(0)
    pdet = chi2_contingency(det.values)[1] if det.shape[1] == 2 and (det.sum(axis=0) > 0).all() else np.nan
    rows.append(dict(family=f, H=H, P=p, H_study=Hs, P_study=ps, P_detection=pdet, n_runs_detected=int((w[f] > 0).sum()),
                     **{f'detected_{h}': (w.loc[g == h, f] > 0).mean() for h in HOSTS}, **{f'median_{h}': w.loc[g == h, f].median() for h in HOSTS},
                     **{f'mean_{h}': w.loc[g == h, f].mean() for h in HOSTS},
                     **{f'study_detected_{h}': (st.loc[st.host_system == h, f] > 0).mean() for h in HOSTS}))
K = pd.DataFrame(rows); K['FDR'] = bh(K.P.values); K['FDR_study'] = bh(K.P_study.values); K['FDR_detection'] = bh(K.P_detection.values)
K.to_csv(f'{O}/metatx_family_tests.csv', index=False); w.to_csv(f'{O}/metatx_rpkm_per_run.csv'); st.to_csv(f'{O}/metatx_rpkm_per_study.csv', index=False)
R['n_sig'] = int((K.FDR < 0.05).sum()); R['sig'] = K.loc[K.FDR < 0.05, 'family'].tolist(); R['n_sig_study'] = int((K.FDR_study < 0.05).sum()); R['sig_study'] = K.loc[K.FDR_study < 0.05, 'family'].tolist()
R['top_host'] = {r.family: max(HOSTS, key=lambda h: getattr(r, f'mean_{h}')) for r in K.itertuples()}
R['targets_expressed'] = int(c[c.NumReads > 0].target.nunique())
json.dump(R, open(f'{O}/metatx_summary.json', 'w'), indent=1, default=float)
print(json.dumps(R, indent=1, default=float))
pd.set_option('display.width', 250)
print(K[['family', 'FDR', 'FDR_study', 'FDR_detection'] + [f'detected_{h}' for h in HOSTS] + [f'median_{h}' for h in HOSTS]].round(4).to_string(index=False))
