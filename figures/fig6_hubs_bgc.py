"""Figure 6. Hub genera and co-location of target genes with biosynthetic gene clusters (BGCs).
usage: python3 fig6_hubs_bgc.py <FINAL dir>
a hub genera (subfamily or family-level labels in >= 3 families) x labelled subfamilies: dot size = curated genes, colour = number
  of metagenome samples with reads on that genus-subfamily pair; b share of genomes with a BGC region overlapping a target gene
  or within 10 kb, hub vs flagship-only genomes; c Fisher odds ratios (95% CI) for each comparison; d co-location per hub genus
  (>= 5 genomes analysed). BGC results: results/s12_bgc (antiSMASH 6.1.1 on the curated set, server).
"""
import sys, os, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import *
import matplotlib.colors as mcolors
B = sys.argv[1]; T = f'{B}/figures/tmp_fig6'; os.makedirs(T, exist_ok=True); out = f'{B}/figures/Figure_6'
t = pd.read_csv(f'{B}/inputs/hitlevel/All_families_taxonomy_curated9396.tsv', sep='\t')
t = t.merge(pd.read_excel(f'{B}/inputs/hitlevel/Domain_architecture_per_hit.xlsx')[['hit_id', 'subfamily_call']], on='hit_id')
t['genus'] = t.classification.fillna('').str.extract(r'g__([^;]*)')[0].fillna('')
lab = t[(t.subfamily_call != 'unassigned') & (t.genus != '')]
nf = lab.groupby('genus').gene.nunique(); hubs = nf[nf >= 3].index
PG = pd.read_csv(f'{B}/results/s12_bgc/per_genus_primary_ge3.tsv', sep='\t')
assert {g.lower() for g in hubs} == {g.lower() for g in PG.genus}, 'hub genera must match the BGC analysis'
FS = pd.read_csv(f'{B}/results/s12_bgc/final_stats.tsv', sep='\t').set_index('set')
GS = pd.read_csv(f'{B}/results/s04_build/genus_subfamily_rpkg_main.csv')

# ---- a: hub genus x subfamily
x = lab[lab.genus.isin(hubs)]
cnt = x.groupby(['genus', 'gene', 'subfamily_call']).size().rename('genes').reset_index()
smp = GS.groupby(['genus', 'family', 'subfamily']).n_samples.max().rename('samples')
cnt['samples'] = [smp.get((g, f, s), 0) for g, f, s in zip(cnt.genus, cnt.gene, cnt.subfamily_call)]
cols = cnt[['gene', 'subfamily_call']].drop_duplicates(); cols['fo'] = cols.gene.map({f: i for i, f in enumerate(FAM)}); cols = cols.sort_values(['fo', 'subfamily_call'])
colkey = list(zip(cols.gene, cols.subfamily_call)); ci = {k: i for i, k in enumerate(colkey)}
rows = cnt.groupby('genus').gene.nunique().sort_values(ascending=False)
rows = sorted(rows.index, key=lambda g: (-rows[g], -cnt[cnt.genus == g].genes.sum())); ri = {g: i for i, g in enumerate(rows)}
NOREAD = '#F4A6A6'
fig = plt.figure(figsize=(174 * MM, 88 * MM)); ax = fig.add_axes([0.14, 0.22, 0.84, 0.54])
norm = mcolors.LogNorm(vmin=1, vmax=max(cnt.samples.max(), 10)); cmap = plt.get_cmap('viridis_r')
for r in cnt.itertuples():
    c = cmap(norm(r.samples)) if r.samples > 0 else NOREAD
    ax.scatter(ci[(r.gene, r.subfamily_call)], ri[r.genus], s=6 + 9 * np.sqrt(r.genes), c=[c], edgecolors='#555555', linewidths=0.3, zorder=3)
ax.set_xlim(-0.6, len(colkey) - 0.4); ax.set_ylim(len(rows) - 0.5, -0.5)
ax.set_xticks(range(len(colkey))); ax.set_xticklabels([s.replace('motifgroup_', 'mg ').replace('no_motif_signal', 'family level') for _, s in colkey], rotation=90, fontsize=5.0)
ax.xaxis.tick_top(); ax.set_yticks(range(len(rows))); ax.set_yticklabels(rows, fontsize=5.2, fontstyle='italic')
for lbl in ax.get_yticklabels():
    if any(ch.isdigit() for ch in lbl.get_text()) or lbl.get_text().isupper(): lbl.set_fontstyle('normal')
ax.grid(color='#EEEEEE', lw=0.4); ax.set_axisbelow(True); ax.tick_params(length=0)
prev = None
for i, (f, s) in enumerate(colkey):
    if f != prev:
        k = sum(1 for ff, _ in colkey if ff == f)
        if i > 0: ax.axvline(i - 0.5, color='#999999', lw=0.5)
        ax.text(i + (k - 1) / 2, -0.015, f, transform=ax.get_xaxis_transform(), ha='center' if k > 2 else 'right', va='top', fontsize=5.4, fontweight='bold', rotation=0 if k > 2 else 60, rotation_mode='anchor')
        prev = f
fig.text(0.14, 0.115, 'Curated genes', fontsize=5.6, fontweight='bold', va='center')
for k_, n_ in enumerate([1, 10, 50]):
    kx = fig.add_axes([0.255 + k_ * 0.055, 0.085, 0.02, 0.06]); kx.axis('off'); kx.set_xlim(0, 1); kx.set_ylim(0, 1)
    kx.scatter([0.5], [0.5], s=6 + 9 * np.sqrt(n_), c='#9ECAE1', edgecolors='#555555', linewidths=0.3, clip_on=False); kx.text(1.5, 0.5, str(n_), va='center', fontsize=5.2)
fig.text(0.50, 0.115, 'Samples with reads', fontsize=5.6, fontweight='bold', va='center')
cax = fig.add_axes([0.66, 0.095, 0.2, 0.025]); cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation='horizontal')
cb.ax.tick_params(labelsize=5, length=2)
nx = fig.add_axes([0.885, 0.02, 0.02, 0.06]); nx.axis('off'); nx.set_xlim(0, 1); nx.set_ylim(0, 1)
nx.scatter([0.5], [0.5], s=24, c=NOREAD, edgecolors='#555555', linewidths=0.3, clip_on=False); nx.text(1.5, 0.5, 'No reads', va='center', fontsize=5.2, fontweight='bold')
fig.savefig(f'{T}/a.pdf'); plt.close(fig)

# ---- b: co-location share
fig, ax = plt.subplots(figsize=(52 * MM, 58 * MM)); fig.subplots_adjust(left=0.24, right=0.98, bottom=0.2, top=0.9)
sets = [('primary_ge3', '≥ 3 families'), ('sensitivity_ge2', '≥ 2 families')]
for k, (sname, lab_) in enumerate(sets):
    r = FS.loc[sname]
    for j, (frac, n, col, c) in enumerate([(r.hub_frac, r.hub_n, r.hub_col, '#0072B2'), (r.fo_frac, r.fo_n, r.fo_col, '#D55E00')]):
        ax.bar(k + (j - 0.5) * 0.36, 100 * frac, width=0.34, color=c)
        ax.text(k + (j - 0.5) * 0.36, 100 * frac + 0.5, f'{int(col)}/{int(n)}', ha='center', va='bottom', fontsize=4.8, rotation=90)
    ax.text(k, 24, pfmt(r.fisher_P) if r.fisher_P < 0.001 else f'P = {r.fisher_P:.2f}', ha='center', fontsize=5.2)
ax.set_xticks(range(2)); ax.set_xticklabels([s[1] for s in sets]); ax.set_ylim(0, 30); ax.set_ylabel('Genomes with co-located BGC (%)')
ax.spines[['top', 'right']].set_visible(False)
from matplotlib.patches import Patch
ax.legend(handles=[Patch(facecolor='#0072B2', label='hub'), Patch(facecolor='#D55E00', label='flagship only')], frameon=False, loc='lower center', bbox_to_anchor=(0.5, 1.0), ncol=2, handlelength=1, fontsize=5.4)
fig.savefig(f'{T}/b.pdf'); plt.close(fig)

# ---- c: odds ratios
order = [('primary_ge3', '≥ 3, 10 kb'), ('primary_ge3_flank5', '≥ 3, 5 kb'), ('primary_ge3_flank20', '≥ 3, 20 kb'), ('primary_ge3_flag4', '≥ 3, original four'),
         ('sensitivity_ge2', '≥ 2, 10 kb'), ('sensitivity_ge2_flag4', '≥ 2, original four')]
fig, ax = plt.subplots(figsize=(60 * MM, 58 * MM)); fig.subplots_adjust(left=0.36, right=0.97, bottom=0.2, top=0.9)
for i, (sname, lab_) in enumerate(order):
    r = FS.loc[sname]; c = '#0072B2' if r.fisher_P < 0.05 else '#D55E00'
    ax.plot([r.fisher_lo, r.fisher_hi], [i, i], color=c, lw=0.9); ax.scatter(r.fisher_OR, i, s=12, c=c, zorder=3)
ax.axvline(1, color='#888888', lw=0.6, ls='--'); ax.set_xscale('log'); ax.set_xticks([0.5, 1, 2, 4, 8]); ax.set_xticklabels(['0.5', '1', '2', '4', '8'])
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter()); ax.set_xlim(0.5, 9)
ax.set_yticks(range(len(order))); ax.set_yticklabels([o[1] for o in order], fontsize=5.4); ax.invert_yaxis()
ax.set_xlabel('Odds ratio (95% CI)'); ax.spines[['top', 'right']].set_visible(False)
fig.savefig(f'{T}/c.pdf'); plt.close(fig)

# ---- d: per genus
P = PG[PG.n_genomes_analysed >= 5].sort_values('fraction', ascending=True)
fig, ax = plt.subplots(figsize=(58 * MM, 58 * MM)); fig.subplots_adjust(left=0.4, right=0.93, bottom=0.2, top=0.97)
ax.barh(range(len(P)), 100 * P.fraction, color='#0072B2', height=0.7)
ax.axvline(100 * FS.loc['primary_ge3', 'fo_frac'], color='#888888', lw=0.6, ls='--')
ax.set_yticks(range(len(P))); ax.set_yticklabels([f'{g} ({n})' for g, n in zip(P.genus, P.n_genomes_analysed)], fontsize=5.2)
for lbl, g in zip(ax.get_yticklabels(), P.genus): lbl.set_fontstyle('normal' if (any(ch.isdigit() for ch in g) or g.isupper()) else 'italic')
ax.set_xlabel('Genomes with co-located BGC (%)'); ax.set_xlim(0, 100); ax.spines[['top', 'right']].set_visible(False)
fig.savefig(f'{T}/d.pdf'); plt.close(fig)

H = assemble([('a', f'{T}/a.pdf', 0, 0, 174, 88), ('b', f'{T}/b.pdf', 0, 91, 52, 58), ('c', f'{T}/c.pdf', 54, 91, 60, 58), ('d', f'{T}/d.pdf', 116, 91, 58, 58)], out)
print('Figure 6 height mm', round(H, 1), '| hubs', len(hubs), '| subfamily columns', len(colkey))
with pd.ExcelWriter(out + '_source_data.xlsx') as xw:
    cnt.to_excel(xw, sheet_name='a_hub_genus_subfamily', index=False); FS.reset_index().to_excel(xw, sheet_name='bc_comparisons', index=False)
    PG.to_excel(xw, sheet_name='d_per_genus', index=False)
