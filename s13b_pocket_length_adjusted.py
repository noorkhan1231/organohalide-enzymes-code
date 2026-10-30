"""Length-adjusted test of top-pocket quality by host system (529 curated structures).
Addresses the reviewer point that pocket counts scale with protein size: pocket metrics are
normalised per 100 residues, and the P2Rank top-3 pocket score is additionally tested on
residuals from a regression on log10 chain length.
usage: python3 s13b_pocket_length_adjusted.py <FINAL dir>"""
import sys, numpy as np, pandas as pd
from scipy.stats import kruskal, spearmanr
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from stats_utils import bh

B = sys.argv[1]
x = pd.ExcelFile(f'{B}/inputs/structures/Supplementary_Data_5.xlsx')
p = x.parse('p2rank_pockets')
keep = set(pd.read_csv(f'{B}/results/s13_structures/structures_529_per_structure.csv').safe_id)
p = p[p.safe_id.isin(keep)].copy()
p['logL'] = np.log10(p['Chain length (residues)'])
MET = {'top3': 'P2Rank mean score, top 3 pockets', 'top1': 'P2Rank score, top pocket',
       'fp100': 'fpocket pockets per 100 residues', 'p2r100': 'P2Rank pockets per 100 residues'}
rows = []
for (fam, sub), d in p.groupby(['Enzyme family', 'Subfamily']):
    if d['Host system'].nunique() < 2 or len(d) < 12:
        continue
    for key, col in MET.items():
        v = d[[col, 'logL', 'Host system']].dropna()
        if v['Host system'].nunique() < 2:
            continue
        groups = [g[col].values for _, g in v.groupby('Host system')]
        H, P = kruskal(*groups)
        # residuals after log10 chain length
        sl, ic = np.polyfit(v.logL, v[col], 1)
        res = v[col] - (sl * v.logL + ic)
        gr = [res[v['Host system'] == h].values for h in sorted(v['Host system'].unique())]
        Hr, Pr = kruskal(*gr)
        rho, prho = spearmanr(v.logL, v[col])
        rows.append(dict(family=fam, subfamily=sub, metric=key, n=len(v),
                         hosts=v['Host system'].nunique(), H=H, P=P,
                         H_length_adjusted=Hr, P_length_adjusted=Pr,
                         rho_length=rho, P_rho=prho))
r = pd.DataFrame(rows)
for m in r.metric.unique():
    i = r.metric == m
    r.loc[i, 'FDR'] = bh(r.loc[i, 'P'].values)
    r.loc[i, 'FDR_length_adjusted'] = bh(r.loc[i, 'P_length_adjusted'].values)
r.round(6).to_csv(f'{B}/results/s13_structures/pocket_length_adjusted_529.csv', index=False)

# tool agreement on the curated set
ag = []
for s, d in p.groupby('Structure set'):
    a = spearmanr(d['fpocket pockets'], d['P2Rank pockets'])
    b = spearmanr(d['fpocket pockets per 100 residues'], d['P2Rank pockets per 100 residues'])
    c = spearmanr(d['Chain length (residues)'], d['fpocket pockets'])
    e = spearmanr(d['Chain length (residues)'], d['P2Rank pockets'])
    ag.append(dict(structure_set=s, n=len(d), rho_counts=a.statistic, P_counts=a.pvalue,
                   rho_per100=b.statistic, P_per100=b.pvalue,
                   rho_length_fpocket=c.statistic, rho_length_p2rank=e.statistic))
a = spearmanr(p['fpocket pockets'], p['P2Rank pockets'])
b = spearmanr(p['fpocket pockets per 100 residues'], p['P2Rank pockets per 100 residues'])
ag.append(dict(structure_set='all 529', n=len(p), rho_counts=a.statistic, P_counts=a.pvalue,
               rho_per100=b.statistic, P_per100=b.pvalue,
               rho_length_fpocket=spearmanr(p['Chain length (residues)'], p['fpocket pockets']).statistic,
               rho_length_p2rank=spearmanr(p['Chain length (residues)'], p['P2Rank pockets']).statistic))
pd.DataFrame(ag).round(4).to_csv(f'{B}/results/s13_structures/pocket_tool_agreement_529.csv', index=False)
print(r[r.metric.isin(['top3', 'top1'])][['family', 'subfamily', 'metric', 'n', 'FDR', 'FDR_length_adjusted', 'rho_length']].to_string(index=False))
print()
print(pd.DataFrame(ag).round(3).to_string(index=False))
