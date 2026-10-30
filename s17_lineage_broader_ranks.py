"""Step 17. Lineage control for the indicators the order-level test could not decide, at broader GTDB ranks.
usage: python3 s17_lineage_broader_ranks.py <FINAL dir> <out dir>
Same procedure as s06_subfamily.py (genome-level presence of the subfamily vs indicator host within strata; CMH odds ratio and
stratified permutation, 4,999 permutations, seed 42), repeated with strata = GTDB family, class and phylum. A broader stratum is a
weaker control (more lineage variation can hide inside it), so a 'host-associated' call at class or phylum is weaker evidence than at order.
BH FDR is applied across all indicators tested at the same rank."""
import sys, os, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stats_utils import bh
from statsmodels.stats.contingency_tables import StratifiedTable
B, out = sys.argv[1:3]; os.makedirs(out, exist_ok=True)
meta = pd.read_csv(f'{B}/inputs/original/target_meta.tsv', sep='\t'); lab = meta[meta.subfamily != 'unassigned'].copy()
sd = pd.read_excel(f'{B}/tables/Supplementary_Data_3.xlsx', sheet_name='genes')[['genome_id', 'GTDB_lineage']].drop_duplicates('genome_id')
rk = {'p': 'phylum', 'c': 'class', 'o': 'order', 'f': 'family'}
for k, v in rk.items(): sd[v + '_x'] = sd.GTDB_lineage.str.extract(rf'{k}__([^;]*)')[0].replace('', np.nan)
lab = lab.merge(sd.rename(columns={'genome_id': 'genome'}), on='genome', how='left')
print('genes without lineage:', int(lab.phylum_x.isna().sum()), 'of', len(lab))
IND = pd.read_excel(f'{B}/figures/Figure_4_source_data.xlsx', sheet_name='c_all_48_indicators')
rng = np.random.default_rng(42)
def test(fname, sname, hosts, col):
    d = lab[lab.family == fname].copy(); d['tax'] = d[col].fillna('NA')
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
R = IND[['family', 'subfamily', 'combination', 'n_genes', 'lineage_order_verdict']].copy()
for lvl, col in [('family', 'family_x'), ('class', 'class_x'), ('phylum', 'phylum_x')]:
    res = [test(r.family, r.subfamily, set(r.combination.split('+')), col) for r in IND.itertuples()]
    R[f'{lvl}_status'] = [x[0] for x in res]; R[f'{lvl}_OR'] = [x[1] for x in res]; R[f'{lvl}_P'] = [x[2] for x in res]; R[f'{lvl}_strata'] = [x[3] for x in res]
    R[f'{lvl}_FDR'] = bh(R[f'{lvl}_P'].values)
    R[f'{lvl}_verdict'] = np.where(R[f'{lvl}_status'] != 'tested', R[f'{lvl}_status'],
                          np.where((R[f'{lvl}_FDR'] < 0.05) & (R[f'{lvl}_OR'] > 1), 'host-associated within lineage', 'not distinguishable from lineage'))
R.to_csv(f'{out}/indicator_lineage_broader_ranks.csv', index=False)
pd.set_option('display.width', 250)
nt = R[R.lineage_order_verdict == 'not testable within lineage']
for lvl in ['family', 'class', 'phylum']: print(lvl, 'among the 19 order-untestable:', nt[f'{lvl}_verdict'].value_counts().to_dict())
for lvl in ['family', 'class', 'phylum']: print(lvl, 'all 48:', R[f'{lvl}_verdict'].value_counts().to_dict())
print(nt[['family', 'subfamily', 'combination', 'n_genes', 'family_verdict', 'class_verdict', 'class_OR', 'class_FDR', 'phylum_verdict', 'phylum_OR', 'phylum_FDR']].round(4).to_string(index=False))
