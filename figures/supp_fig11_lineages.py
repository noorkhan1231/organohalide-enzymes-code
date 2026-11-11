"""Supplementary Figure 16 (added 6 Nov 2026, Prof. Zhang comment C38): bacterial orders carrying the six host-linked subfamilies.
usage: python3 supp_fig11_lineages.py <FINAL dir>
Rows are the GTDB orders of all genomes that carry one of the six subfamilies that remain associated with their host system
within GTDB order (lineage control, Table S12), arranged as a GTDB taxonomy cladogram (phylum > class > order; branch lengths
are not evolutionary distances). Columns are the six subfamilies; in each cell, one circle per host system, area = carrier genomes.
"""
import sys, os, re, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from figstyle import *
B = sys.argv[1]; OUT = f'{B}/figures/supplementary'
g = pd.read_excel(f'{B}/tables/Supplementary_Data_3.xlsx', sheet_name='genes')
SIX = [('GSTDase', 'motifgroup_SxC', 'GSTDase mg SxC', 'Sponge'), ('GSTDase', 'linD', 'GSTDase linD', 'Sponge'),
       ('HyDase', 'dhmA', 'HyDase dhmA', 'Rhizosphere'), ('HyDase', 'motifgroup_DRG', 'HyDase mg DRG', 'Chicken'),
       ('VanHase', 'V-HPO', 'VanHase V-HPO', 'Coral + Rhizosphere'), ('VanHase', 'V-IPO_ambiguous', 'VanHase V-IPO_amb.', 'Coral + Rhizosphere')]
key = {(f, s): lab for f, s, lab, _ in SIX}
g = g[[(f, s) in key for f, s in zip(g.enzyme_family, g.subfamily)]].copy()
g['lab'] = [key[(f, s)] for f, s in zip(g.enzyme_family, g.subfamily)]
rk = lambda k: g.GTDB_lineage.str.extract(rf'{k}__([^;]*)')[0].replace('', np.nan).fillna('unclassified')
g['phylum'], g['class'], g['order'] = rk('p'), rk('c'), rk('o')
G = g.drop_duplicates(['genome_id', 'lab'])
cnt = G.groupby(['phylum', 'class', 'order', 'lab', 'host_system']).genome_id.nunique().reset_index(name='genomes')
# order of rows: phyla by number of genomes, classes and orders likewise
tot = G.groupby(['phylum', 'class', 'order']).genome_id.nunique().reset_index(name='n')
pn = tot.groupby('phylum').n.sum(); cn = tot.groupby(['phylum', 'class']).n.sum()
tot['kp'] = tot.phylum.map(-pn); tot['kc'] = [-cn[(p, c)] for p, c in zip(tot.phylum, tot['class'])]
tot = tot.sort_values(['kp', 'phylum', 'kc', 'class', 'n', 'order'], ascending=[True, True, True, True, False, True]).reset_index(drop=True)
rows = list(zip(tot.phylum, tot['class'], tot.order)); nR = len(rows)
H = 14 + nR * 2.6
fig = plt.figure(figsize=(174 * MM, H * MM))
ax = fig.add_axes([0.005, 6 / H, 0.99, (H - 8) / H]); ax.set_xlim(0, 174); ax.set_ylim(nR - 0.5, -5.2); ax.axis('off')
# cladogram: x positions of levels (mm)
XR, XP, XC, XO, XL = 2, 30, 41, 51, 52.5          # root, phylum, class, order nodes, label start
def vline(x, y0, y1): ax.plot([x, x], [y0, y1], color='#555555', lw=0.45)
def hline(x0, x1, y): ax.plot([x0, x1], [y, y], color='#555555', lw=0.45)
ys = {r: i for i, r in enumerate(rows)}
py = {}
for p in dict.fromkeys(tot.phylum):
    cy = {}
    for c in dict.fromkeys(tot[tot.phylum == p]['class']):
        oy = [ys[r] for r in rows if r[0] == p and r[1] == c]
        for y in oy: hline(XC, XO, y)
        vline(XC, min(oy), max(oy)); cy[c] = (min(oy) + max(oy)) / 2; hline(XP, XC, cy[c])
    vline(XP, min(cy.values()), max(cy.values())); py[p] = (min(cy.values()) + max(cy.values())) / 2; hline(XR, XP, py[p])
    ax.text(XP - 0.6, py[p] - 0.25, p, fontsize=4.3, color='#333333', va='bottom', ha='right', style='italic')
vline(XR, min(py.values()), max(py.values()))
for (p, c, o), y in ys.items():
    ax.text(XL, y, o, fontsize=4.6, va='center', ha='left')
# matrix
X0, DX = 91, 13.6
for j, (f, s, lab, ind) in enumerate(SIX):
    x = X0 + j * DX
    sub = G[G.lab == lab]
    nm = lab.split(' ', 1)
    ax.text(x + 4.5, -5.0, nm[0], fontsize=5.2, ha='center', va='bottom', fontweight='bold')
    ax.text(x + 4.5, -4.1, nm[1], fontsize=5.2, ha='center', va='bottom', fontweight='bold')
    ax.text(x + 4.5, -2.75, f'{sub.genome_id.nunique()} genomes', fontsize=4.2, ha='center', va='bottom', color='#444444')
    ax.text(x + 4.5, -1.95, f'{sub.order.nunique()} orders, {sub.phylum.nunique()} phyla', fontsize=4.2, ha='center', va='bottom', color='#444444')
    ax.text(x + 4.5, -1.15, ind.replace('Coral + Rhizosphere', 'Co + Rh'), fontsize=4.2, ha='center', va='bottom', color='#444444')
    for k, h in enumerate(HOSTS): ax.text(x + k * 3.0, -0.55, h[:2], fontsize=3.6, ha='center', va='bottom', color=COL[h])
    for y in range(nR):
        if y % 2 == 0: ax.add_patch(matplotlib.patches.Rectangle((x - 1.4, y - 0.5), 11.8, 1, color='#F4F4F4', lw=0, zorder=0))
    for k, h in enumerate(HOSTS):
        d = cnt[(cnt.lab == lab) & (cnt.host_system == h)]
        for r in d.itertuples():
            ax.scatter(x + k * 3.0, ys[(r.phylum, r[2], r.order)], s=3.2 * np.sqrt(r.genomes) ** 2 if False else 2.2 + 5.5 * np.sqrt(r.genomes), color=COL[h], lw=0, zorder=3)
for y in range(nR): ax.plot([XL + 23, X0 - 1.8], [y, y], color='#DDDDDD', lw=0.25, zorder=0)
# legend
lx = 6
for k, h in enumerate(HOSTS):
    fig.text((lx + k * 22) / 174, 2.4 / H, '●', color=COL[h], fontsize=6, va='center'); fig.text((lx + 3 + k * 22) / 174, 2.4 / H, h, fontsize=5.2, va='center')
for k, n in enumerate([1, 10, 50]):
    xx = (100 + k * 16) / 174
    fig.add_artist(matplotlib.lines.Line2D([xx], [2.4 / H], marker='o', color='#777777', markersize=np.sqrt(2.2 + 5.5 * np.sqrt(n)) * 0.95, lw=0, transform=fig.transFigure))
    fig.text(xx + 3 / 174, 2.4 / H, f'{n} genome' + ('s' if n > 1 else ''), fontsize=5.0, va='center')
fig.text(0.012, (H - 3) / H, 'GTDB phylum › class › order', fontsize=5.2, fontweight='bold', va='center')
save(fig, f'{OUT}/Supplementary_Figure_11'); plt.close(fig)
with pd.ExcelWriter(f'{OUT}/Supplementary_Figure_11_source_data.xlsx') as xw:
    cnt.rename(columns={'lab': 'subfamily'}).to_excel(xw, sheet_name='carrier_genomes', index=False)
    G.groupby('lab').agg(genomes=('genome_id', 'nunique'), orders=('order', 'nunique'), classes=('class', 'nunique'), phyla=('phylum', 'nunique')).reset_index().rename(columns={'lab': 'subfamily'}).to_excel(xw, sheet_name='summary', index=False)
print('SF16 rows (orders):', nR, 'phyla:', tot.phylum.nunique(), 'height mm', round(H, 1))
