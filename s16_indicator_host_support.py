"""Step 16. Host support for indicator subfamilies: surroundings contrast and expression, at subfamily level.
usage: python3 s16_indicator_host_support.py <inputs_dir> <results_dir> <out_dir>

For each of the 48 indicator subfamilies (Figure 4 source data, sheet c_all_48_indicators), and in the indicator host system(s):
 (1) Host vs surroundings: the s11 procedure applied to the subfamily instead of the family. Host RPKG at sample level (s04
     subfamily_rpkg_long_main), surroundings RPKG at run level; null-normalised abundance on runs with a null value (s11
     ambient_runs_rpkg null_RPKG). Two-sided Mann-Whitney U on study means, >= 3 studies per side; BH FDR across all
     subfamily x system tests within each measure. Call as in s11: host-enriched only when RPKG and null-normalised agree
     in direction and both reach FDR < 0.05.
 (2) Expression: subfamily RPKM in the 190 metatranscriptome runs of s07; share of runs and of studies with reads per host
     system; Kruskal-Wallis across host systems at run level and at study level (BH FDR across subfamilies).
Nothing upstream is changed; this step only reads earlier outputs and the per-target counts used by s11 and s07.
"""
import sys, os, csv, json, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stats_utils import HOSTS, bh, kw
from scipy.stats import mannwhitneyu
I, RES, out = sys.argv[1:4]; os.makedirs(out, exist_ok=True)
FIG = os.path.join(os.path.dirname(RES), 'figures', 'Figure_4_source_data.xlsx')
IND = pd.read_excel(FIG, sheet_name='c_all_48_indicators')
meta = pd.read_csv(f'{I}/original/target_meta.tsv', sep='\t')[['target', 'family', 'subfamily']]

# ---------------- (1) surroundings
R = pd.read_csv(f'{RES}/s11_ambient/ambient_runs_rpkg.csv', dtype={'run': str})
cnt = [pd.read_csv(f'{I}/original/per_target_counts.tsv.gz', sep='\t', dtype={'run': str})]
for j in ['job6', 'job6b', 'job14', 'job16']:
    cnt.append(pd.read_csv(f'{I}/{j}/extra_counts.tsv.gz', sep='\t', dtype={'run': str}))
c = pd.concat(cnt, ignore_index=True); c = c[c.dataset != 'metatx']
c['ds'] = np.where(c.dataset == 'ambient', 'ambient', 'host')
c = c.merge(meta, on=['target', 'family'], how='inner')
c['RPK'] = c.NumReads / (c.Length / 1000)
key = c.family + '|' + c.subfamily
c['key'] = key
sub_run = c.groupby(['ds', 'run', 'key']).RPK.sum()
hs = pd.read_csv(f'{RES}/s04_build/subfamily_rpkg_long_main.csv')
hs['key'] = hs.family + '|' + hs.subfamily
fam_main = pd.read_csv(f'{RES}/s04_build/family_rpkg_main.csv')[['sample_id', 'host_system', 'study']]

def proj(h, a):
    d = dict(n_host_studies=len(h), n_ambient_studies=len(a))
    if len(h) < 3 or len(a) < 3: return d
    p = mannwhitneyu(h, a, alternative='two-sided').pvalue
    return dict(d, mean_host=h.mean(), mean_ambient=a.mean(),
                ratio=h.mean() / a.mean() if a.mean() > 0 else (np.inf if h.mean() > 0 else np.nan),
                host_studies_detected=int((h > 0).sum()), ambient_studies_detected=int((a > 0).sum()), P=p)

T = []
for r in IND.itertuples():
    k = f'{r.family}|{r.subfamily}'
    hsk = hs[hs.key == k].set_index('sample_id').RPKG
    for sysn in r.combination.split('+'):
        # RPKG: host samples (all main-set samples of the system, zero where absent) vs surrounding runs
        H = fam_main[fam_main.host_system == sysn].copy(); H['v'] = H.sample_id.map(hsk).fillna(0)
        A = R[(R.system == sysn) & (R.group != 'host')].copy()
        A['v'] = [sub_run.get(('ambient' if d == 'ambient' else 'host', rn, k), 0.0) / ge for d, rn, ge in zip(A.ds, A.run, A.GE)]
        t1 = proj(H.groupby('study').v.mean(), A.groupby('unit').v.mean())
        # null-normalised: runs with a null value
        X = R[(R.system == sysn) & (R.null_RPKG > 0)].copy()
        X['v'] = [sub_run.get((d, rn, k), 0.0) / ge for d, rn, ge in zip(X.ds, X.run, X.GE)]
        X['norm'] = X.v / X.null_RPKG
        t2 = proj(X[X.group == 'host'].groupby('unit').norm.mean(), X[X.group != 'host'].groupby('unit').norm.mean())
        T.append(dict(family=r.family, subfamily=r.subfamily, combination=r.combination, system=sysn, n_genes=r.n_genes,
                      lineage_order=r.lineage_order_verdict,
                      **{f'{a}_RPKG': b for a, b in t1.items()}, **{f'{a}_norm': b for a, b in t2.items()}))
T = pd.DataFrame(T)
T['FDR_RPKG'] = np.nan; T['FDR_norm'] = np.nan
m = T.P_RPKG.notna(); T.loc[m, 'FDR_RPKG'] = bh(T.loc[m, 'P_RPKG'].values)
m = T.P_norm.notna(); T.loc[m, 'FDR_norm'] = bh(T.loc[m, 'P_norm'].values)
def call(x):
    if not (x.FDR_RPKG < 0.05 and x.FDR_norm < 0.05): return 'no consistent difference' if pd.notna(x.P_RPKG) else 'not testable'
    if x.ratio_RPKG > 1 and x.ratio_norm > 1: return 'host-enriched'
    if x.ratio_RPKG < 1 and x.ratio_norm < 1: return 'surroundings-enriched'
    return 'unresolved'
T['call'] = T.apply(call, axis=1)
T.to_csv(f'{out}/indicator_host_vs_surroundings.csv', index=False)

# ---------------- (2) expression
W = pd.read_csv(f'{RES}/s07_metatx/metatx_rpkm_per_run.csv')[['run_accession', 'host_system', 'study', 'million']]
mc = pd.read_csv(f'{I}/metatx/metatx_counts.tsv.gz', sep='\t'); mc = mc[mc.dataset == 'metatx']
mc = mc[mc.run.isin(W.run_accession)].merge(meta, on=['target', 'family'], how='inner')
mc['key'] = mc.family + '|' + mc.subfamily; mc['rpk'] = mc.NumReads / (mc.Length / 1000)
sm = mc.groupby(['run', 'key']).rpk.sum()
E = []
for r in IND.itertuples():
    k = f'{r.family}|{r.subfamily}'
    w = W.copy(); w['v'] = [sm.get((rn, k), 0.0) / mi for rn, mi in zip(w.run_accession, w.million)]
    st = w.groupby(['host_system', 'study']).v.mean().reset_index()
    H, p = kw(w.v.values, w.host_system.values) if (w.v > 0).any() else (np.nan, np.nan)
    Hs, ps = kw(st.v.values, st.host_system.values) if (st.v > 0).any() else (np.nan, np.nan)
    d = dict(family=r.family, subfamily=r.subfamily, combination=r.combination, n_genes=r.n_genes, lineage_order=r.lineage_order_verdict,
             runs_detected=int((w.v > 0).sum()), P_run=p, P_study=ps)
    for h in HOSTS:
        d[f'runs_{h}'] = int((w.host_system == h).sum()); d[f'det_runs_{h}'] = float((w.loc[w.host_system == h, 'v'] > 0).mean())
        d[f'det_studies_{h}'] = int((st.loc[st.host_system == h, 'v'] > 0).sum()); d[f'studies_{h}'] = int((st.host_system == h).sum())
        d[f'mean_{h}'] = float(w.loc[w.host_system == h, 'v'].mean())
    d['top_host'] = max(HOSTS, key=lambda h: d[f'mean_{h}']) if d['runs_detected'] else ''
    d['top_in_indicator'] = d['top_host'] in r.combination.split('+')
    E.append(d)
E = pd.DataFrame(E)
for col in ['P_run', 'P_study']:
    E['FDR' + col[1:]] = np.nan; m = E[col].notna(); E.loc[m, 'FDR' + col[1:]] = bh(E.loc[m, col].values)
E.to_csv(f'{out}/indicator_expression.csv', index=False)

pd.set_option('display.width', 250); pd.set_option('display.max_columns', 40)
six = T[T.lineage_order == 'host-associated within lineage']
print(six[['family', 'subfamily', 'system', 'n_host_studies_RPKG', 'n_ambient_studies_RPKG', 'host_studies_detected_RPKG', 'ambient_studies_detected_RPKG',
           'ratio_RPKG', 'FDR_RPKG', 'ratio_norm', 'FDR_norm', 'call']].round(4).to_string(index=False))
print(T.groupby(['lineage_order', 'call']).size().unstack(fill_value=0))
e6 = E[E.lineage_order == 'host-associated within lineage']
print(e6[['family', 'subfamily', 'combination', 'runs_detected', 'FDR_run', 'FDR_study', 'top_host', 'top_in_indicator'] +
         [f'det_runs_{h}' for h in HOSTS] + [f'det_studies_{h}' for h in HOSTS]].round(3).to_string(index=False))
print(E.groupby('lineage_order').agg(n=('subfamily', 'size'), expressed=('runs_detected', lambda v: int((v > 0).sum())),
                                     top_in_ind=('top_in_indicator', 'sum'), sig_run=('FDR_run', lambda v: int((v < 0.05).sum()))))
json.dump(dict(surroundings=T.groupby('call').size().to_dict(), expression_sig_run=int((E.FDR_run < 0.05).sum()),
               expression_sig_study=int((E.FDR_study < 0.05).sum()), expressed=int((E.runs_detected > 0).sum())),
          open(f'{out}/summary.json', 'w'), indent=1)
