"""Step 9. Subfamily-level statistics on sample-level RPKG (all targets; every community sample kept, E4 fix).
usage: python3 s06_subfamily.py <subfamily_rpkg_long_<set>.csv> <family_rpkg_<set>.csv> <target_meta.tsv> <outdir> [label]
Per family: PERMANOVA/PERMDISP on subfamily composition (sample and study level), Kruskal-Wallis per subfamily,
indicator subfamilies (multipatt, 999 permutations, BH FDR within family), genome-level lineage control (CMH within GTDB
order or genus strata, stratified permutation), subfamily co-occurrence."""
import sys, os, json, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lib'))
from stats_utils import *
from multipatt import multipatt
from statsmodels.stats.contingency_tables import StratifiedTable
sub_f, fam_f, meta_f, out = sys.argv[1:5]; label = sys.argv[5] if len(sys.argv) > 5 else 'all'
os.makedirs(out, exist_ok=True); P = f'{out}/{label}'
fam = pd.read_csv(fam_f); fam = fam[fam[FAM].sum(axis=1) > 0].reset_index(drop=True)      # community samples
s = pd.read_csv(sub_f); s = s[s.sample_id.isin(fam.sample_id) & (s.subfamily != 'unassigned')]
meta = pd.read_csv(meta_f, sep='\t')
RES = {}; KW, IND, PERM = [], [], []
SUBFAMS = [f for f in ['RDase', 'HyDase', 'GSTDase', 'OxDase', 'MetDase', 'FlaHase', 'VanHase', 'NHFeHase'] if f in s.family.unique()]
mats = {}
for f in SUBFAMS:
    w = s[s.family == f].pivot_table(index='sample_id', columns='subfamily', values='RPKG', aggfunc='sum', fill_value=0)
    w = fam[['sample_id', 'host_system', 'study']].set_index('sample_id').join(w).fillna(0)
    cols = [c for c in w.columns if c not in ('host_system', 'study')]; mats[f] = (w, cols)
    g = w.host_system.values
    # composition: PERMANOVA/PERMDISP on samples with any labelled reads (Bray-Curtis undefined otherwise)
    if len(cols) >= 2:
        V = run_vegan(w.reset_index()[['sample_id', 'host_system', 'study'] + cols], 'sample_id', 'host_system', 'study', '-', f'{P}_{f}_vegan', int(os.environ.get('NMDS_K', '0')))
        v = V.set_index('test')
        PERM.append(dict(family=f, n_subfamilies=len(cols), n_samples=int(v.loc['host_sample', 'n']), R2=v.loc['host_sample', 'R2'], P=v.loc['host_sample', 'P'],
                         permdisp_F=v.loc['permdisp_sample', 'F'], permdisp_P=v.loc['permdisp_sample', 'P'],
                         R2_study=v.loc['host_study', 'R2'] if 'host_study' in v.index else np.nan, P_study=v.loc['host_study', 'P'] if 'host_study' in v.index else np.nan,
                         nmds_stress=pd.read_csv(f'{P}_{f}_vegan_nmds.csv').stress.iat[0] if os.path.exists(f'{P}_{f}_vegan_nmds.csv') else np.nan))
    # abundance of each subfamily among all community samples
    for c in cols:
        x = w[c].values; npos = int((x > 0).sum())
        H, p = kw(x, g) if npos >= 10 else (np.nan, np.nan)
        KW.append(dict(family=f, subfamily=c, n_detected=npos, H=H, P=p, **{f'mean_{h}': x[g == h].mean() for h in HOSTS}, **{f'prev_{h}': (x[g == h] > 0).mean() for h in HOSTS}))
    # indicators among all community samples
    lev, comb, stat, pv = multipatt(w[cols].values, list(g), 999, 42)
    r = pd.DataFrame(dict(family=f, subfamily=cols, combination=['+'.join(c) for c in comb], stat=stat, P=pv)); r['FDR'] = bh(r.P.values)
    for h in HOSTS: r[f'prev_{h}'] = [(w.loc[g == h, c] > 0).mean() for c in cols]
    IND.append(r)
K = pd.DataFrame(KW); K['FDR'] = np.nan
for f in K.family.unique():
    k = K.family == f; K.loc[k, 'FDR'] = bh(K.loc[k, 'P'].values)
K.to_csv(f'{P}_subfamily_kruskal.csv', index=False)
I = pd.concat(IND); I.to_csv(f'{P}_subfamily_indicators_all.csv', index=False)
Pm = pd.DataFrame(PERM); Pm.to_csv(f'{P}_subfamily_permanova.csv', index=False)
sig = I[I.FDR <= 0.05].copy()
# genes and genomes behind each indicator
lab = meta[meta.subfamily != 'unassigned']
gc = lab.groupby(['family', 'subfamily']).agg(n_genes=('target', 'size'), n_genomes=('genome', 'nunique')).reset_index()
hc = lab.groupby(['family', 'subfamily', 'origin_host']).size().unstack(fill_value=0).reindex(columns=HOSTS, fill_value=0).add_prefix('genes_').reset_index()
sig = sig.merge(gc, on=['family', 'subfamily'], how='left').merge(hc, on=['family', 'subfamily'], how='left')
sig['single_host'] = ~sig.combination.str.contains(r'\+')

# lineage control (genome level): among genomes carrying labelled genes of the family, is carriage of the indicator
# subfamily associated with the indicated host(s) within GTDB order (or family) strata?  CMH + stratified permutation.
rng = np.random.default_rng(42)
def lineage_test(fname, sname, hosts, level):
    d = lab[lab.family == fname].copy()
    d['tax'] = d[level].fillna('NA')
    gg = d.groupby('genome').agg(carry=('subfamily', lambda x: int((x == sname).any())), host=('origin_host', 'first'), stratum=('tax', 'first')).reset_index()
    gg['inH'] = gg.host.isin(hosts).astype(int)
    if gg.carry.sum() < 3: return 'too few carrier genomes', np.nan, np.nan, 0
    st = [k for k, x in gg.groupby('stratum') if x.carry.nunique() == 2 and x.inH.nunique() == 2]
    ok = gg[gg.stratum.isin(st)]
    if len(st) < 1 or ok.carry.sum() < 3: return 'not testable within lineage', np.nan, np.nan, len(st)
    tabs = [np.array([[((x.carry == 1) & (x.inH == 1)).sum(), ((x.carry == 1) & (x.inH == 0)).sum()],
                      [((x.carry == 0) & (x.inH == 1)).sum(), ((x.carry == 0) & (x.inH == 0)).sum()]]) for _, x in ok.groupby('stratum')]
    try: OR = StratifiedTable(tabs).oddsratio_pooled
    except Exception: OR = np.nan
    obs = sum(t[0, 0] for t in tabs); grp = [(x.carry.values, x.inH.values) for _, x in ok.groupby('stratum')]
    cnt = sum(sum((c * rng.permutation(i)).sum() for c, i in grp) >= obs for _ in range(4999))
    return 'tested', OR, (cnt + 1) / 5000, len(st)
for lvl in ['order', 'genus']:
    res = [lineage_test(r.family, r.subfamily, set(r.combination.split('+')), lvl) for r in sig.itertuples()]
    sig[f'lineage_{lvl}_status'] = [x[0] for x in res]; sig[f'lineage_{lvl}_OR'] = [x[1] for x in res]; sig[f'lineage_{lvl}_P'] = [x[2] for x in res]; sig[f'lineage_{lvl}_strata'] = [x[3] for x in res]
    sig[f'lineage_{lvl}_FDR'] = bh(sig[f'lineage_{lvl}_P'].values)
    sig[f'lineage_{lvl}_verdict'] = np.where(sig[f'lineage_{lvl}_status'] != 'tested', sig[f'lineage_{lvl}_status'],
                                    np.where((sig[f'lineage_{lvl}_FDR'] < 0.05) & (sig[f'lineage_{lvl}_OR'] > 1), 'host-associated within lineage', 'not distinguishable from lineage'))
sig.to_csv(f'{P}_subfamily_indicators_significant.csv', index=False)
# subfamily co-occurrence network (all community samples; subfamilies detected in >= 10 samples)
allw = pd.concat([mats[f][0][mats[f][1]].add_prefix(f + '::') for f in SUBFAMS], axis=1)
allw = allw.loc[:, (allw > 0).sum() >= 10]
C = cooccur(allw.values, list(allw.columns)); C.to_csv(f'{P}_subfamily_cooccurrence.csv', index=False)
Cs = C[(C.fdr < 0.05) & (C.rho.abs() > 0.3)]
FAMIND = {'RDase', 'OxDase', 'DehyDase', 'VanHase', 'NHFeHase', 'SAMHase'}
RES = dict(n_indicators=int(len(sig)), by_family=sig.family.value_counts().to_dict(), single_host=int(sig.single_host.sum()),
           indicators_with_le3_genes=int((sig.n_genes <= 3).sum()), cooc_pairs=int(len(Cs)), cooc_negative=int((Cs.rho < 0).sum()),
           lineage_order=sig.lineage_order_verdict.value_counts().to_dict(), lineage_genus=sig.lineage_genus_verdict.value_counts().to_dict(),
           permanova=Pm[['family', 'n_samples', 'R2', 'P', 'permdisp_P', 'R2_study', 'P_study']].round(4).values.tolist())
json.dump(RES, open(f'{P}_subfamily_summary.json', 'w'), indent=1, default=float)
print(json.dumps(RES, indent=1, default=float))
pd.set_option('display.width', 250)
print(sig[['family', 'subfamily', 'combination', 'stat', 'FDR', 'n_genes', 'n_genomes'] + [f'prev_{h}' for h in HOSTS] + ['lineage_order_verdict']].round(3).to_string(index=False))
