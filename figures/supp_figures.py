"""Supplementary Figures 4-15 (new data) and copies of Supplementary Figures 1-3 (gene-level, unchanged).
usage: python3 supp_figures.py <FINAL dir> <old package dir>
SF1 E-value cut-offs, SF2 gene trees, SF3 prevalence per carrier genome: copied from the earlier package (no sample dependence).
SF4 dataset (genome equivalents, reads, samples per country); SF5 robustness (genome-equivalent thresholds, leave-one-study-out);
SF6 family abundance at study level; SF7 subfamily ordinations (PCoA, Bray-Curtis); SF8 subfamily co-occurrence;
SF9 geography (Mantel); SF10 host vs surroundings per family; SF11 metatranscriptomes at study level; SF12 signal peptides;
SF13 structural PERMANOVA and PCoA; SF14 chain length and pockets of four subfamilies; SF15 breadth of enzyme repertoire per carrier genus.
"""
import sys, os, re, shutil, json, glob, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from figstyle import *
from matplotlib.patches import Patch
from matplotlib.legend_handler import HandlerTuple
def tup(hatch=None): return tuple(Patch(facecolor=('white' if hatch else COL[h]), edgecolor=COL[h], hatch=hatch, lw=0.6) for h in HOSTS)
from scipy.spatial.distance import pdist, squareform
from stats_utils import kw
B, OLD = sys.argv[1:3]; R = f'{B}/results'; OUT = f'{B}/figures/supplementary'; os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(21)
def pcoa(D):
    n = D.shape[0]; J = np.eye(n) - 1 / n; Bm = -0.5 * J @ (D ** 2) @ J; w, v = np.linalg.eigh(Bm); i = np.argsort(w)[::-1]; w, v = w[i], v[:, i]
    pos = w > 0; return v[:, :2] * np.sqrt(np.maximum(w[:2], 0)), w[:2] / w[pos].sum()
def logticks(ax, lo, hi, axis='y'):
    t = range(lo, hi + 1); lab = [f'$10^{{{k}}}$' for k in t]
    (ax.set_yticks if axis == 'y' else ax.set_xticks)(list(t)); (ax.set_yticklabels if axis == 'y' else ax.set_xticklabels)(lab)
def strip_box(ax, vals, pos, col, w=0.5, s=1.5):
    if len(vals) >= 3:
        ax.boxplot(vals, positions=[pos], widths=w, showfliers=False, patch_artist=True, medianprops=dict(color='black', lw=0.6),
                   boxprops=dict(facecolor=col, alpha=0.55, lw=0.4), whiskerprops=dict(lw=0.4), capprops=dict(lw=0.4))
    ax.scatter(pos + rng.uniform(-w / 3, w / 3, len(vals)), vals, s=s, c=col, linewidths=0, zorder=3, rasterized=True)

# ---- SF1-3 copies
import pymupdf
for k in [1, 2, 3]:   # gene-level figures from the earlier package; panel letters converted to uppercase (Cell Press style)
    d = pymupdf.open(f'{OLD}/Supplementary_figures/Supplementary_Figure_{k}.pdf'); pg = d[0]; spans = []
    for b in pg.get_text('dict')['blocks']:
        for l in b.get('lines', []):
            for sp in l['spans']:
                t = sp['text'].strip()
                if len(t) == 1 and t.isalpha() and t.islower() and 'Bold' in sp['font'] and sp['size'] >= 6.5: spans.append((t, sp['bbox'], sp['origin']))
    for t, bb, org in spans: pg.add_redact_annot(pymupdf.Rect(bb), fill=(1, 1, 1))
    pg.apply_redactions(images=pymupdf.PDF_REDACT_IMAGE_NONE)
    for t, bb, org in spans: pg.insert_text((org[0] - 1.2, org[1]), t.upper(), fontsize=7.5, fontname='LSB', fontfile='/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf')
    d.save(f'{OUT}/Supplementary_Figure_{k}.pdf', garbage=3, deflate=True)
    pg.get_pixmap(dpi=1200).save(f'{OUT}/Supplementary_Figure_{k}.jpg', jpg_quality=92); pg.get_pixmap(dpi=200).save(f'{OUT}/Supplementary_Figure_{k}_preview.png')
    print(f'Supplementary Figure {k}: {len(spans)} panel letters converted')

S = pd.read_csv(f'{R}/s03_qc/samples.tsv', sep='\t'); S = S[S.set_main]
F = pd.read_csv(f'{R}/s04_build/family_rpkg_main.csv')

# ---- SF4 dataset
CO = pd.read_csv(f'{R}/s08_map/sample_coordinates.tsv', sep='\t')
fig = plt.figure(figsize=(174 * MM, 120 * MM))
ax = fig.add_axes([0.07, 0.6, 0.38, 0.36])
for j, h in enumerate(HOSTS): strip_box(ax, np.log10(S[S.host_system == h].GE.clip(lower=0.1)), j, COL[h], 0.6, 1.2)
ax.set_xticks(range(4)); ax.set_xticklabels([f'{h}\n({(S.host_system == h).sum()})' for h in HOSTS]); logticks(ax, -1, 5); ax.set_ylabel('Genome equivalents per sample')
ax.spines[['top', 'right']].set_visible(False); label(ax, 'a', -0.15)
ax = fig.add_axes([0.58, 0.6, 0.38, 0.36])
for j, h in enumerate(HOSTS): strip_box(ax, np.log10(S[S.host_system == h].reads_processed), j, COL[h], 0.6, 1.2)
ax.set_xticks(range(4)); ax.set_xticklabels([f'{h}\n({(S.host_system == h).sum()})' for h in HOSTS]); logticks(ax, 5, 9); ax.set_ylabel('Reads processed per sample')
ax.spines[['top', 'right']].set_visible(False); label(ax, 'b', -0.15)
c = CO.country.fillna('not recorded').replace({'United States of America': 'USA', 'Viet Nam': 'Vietnam', 'United Arab Emirates': 'UAE', 'United Kingdom': 'UK'})
Tc = pd.crosstab(c, CO.host_system).reindex(columns=HOSTS, fill_value=0); Tc = Tc.loc[Tc.sum(axis=1).sort_values(ascending=False).index]
nr = int(Tc.loc['not recorded'].sum()) if 'not recorded' in Tc.index else 0; Tc = Tc[Tc.index != 'not recorded']
ax = fig.add_axes([0.07, 0.13, 0.91, 0.33]); bottom = np.zeros(len(Tc))
for h in HOSTS: ax.bar(range(len(Tc)), Tc[h].values, bottom=bottom, color=COL[h], width=0.75, linewidth=0, label=h); bottom += Tc[h].values
ax.set_xticks(range(len(Tc))); ax.set_xticklabels(Tc.index, rotation=60, ha='right', rotation_mode='anchor', fontsize=4.6); ax.set_xlim(-0.6, len(Tc) - 0.4)
ax.set_ylabel('Samples'); ax.spines[['top', 'right']].set_visible(False); ax.legend(ncol=4, frameon=False, loc='upper right', handlelength=1); label(ax, 'c', -0.06)
save(fig, f'{OUT}/Supplementary_Figure_4'); plt.close(fig)
with pd.ExcelWriter(f'{OUT}/Supplementary_Figure_4_source_data.xlsx') as xw:
    S[['sample_id', 'host_system', 'GE', 'reads_processed']].to_excel(xw, sheet_name='ab', index=False); Tc.to_excel(xw, sheet_name='c_samples_per_country')

# ---- SF5 robustness
SE = pd.read_csv(f'{R}/s05_family_sens/sensitivity_summary.csv'); LOO = pd.read_csv(f'{R}/s05_family/main_loo_kw.csv')
fig = plt.figure(figsize=(174 * MM, 62 * MM))
ax = fig.add_axes([0.07, 0.2, 0.25, 0.68]); lab = {'main': 'all', 'ge20': '≥ 20', 'ge50': '≥ 50', 'ge200': '≥ 200'}
x = np.arange(len(SE)); ax.bar(x - 0.18, SE.R2, 0.34, color='#3B5B8C', label='sample level'); ax.bar(x + 0.18, SE.R2_study, 0.34, color='#9AAFCB', label='study level')
ax.set_xticks(x); ax.set_xticklabels([f'{lab[s]}\n{n:,}' for s, n in zip(SE.set, SE.n)], fontsize=5.4); ax.set_xlabel('Genome equivalents (samples)', fontsize=6)
for xi, (p1, p2, r1, r2_) in enumerate(zip(SE.P, SE.P_study, SE.R2, SE.R2_study)):
    ax.text(xi - 0.18, r1 + 0.006, pfmt(p1).replace('P < 0.001', 'P ≤ 0.001') if p1 > 0.001 else 'P = 0.001', rotation=90, ha='center', va='bottom', fontsize=4.6)
    ax.text(xi + 0.18, r2_ + 0.006, 'P = 0.0001' if p2 <= 0.0001 else f'P = {p2:.3f}', rotation=90, ha='center', va='bottom', fontsize=4.6)
ax.set_ylabel('R² (host system)'); ax.set_ylim(0, 0.38); ax.legend(frameon=False, loc='upper left', ncol=2, fontsize=5.4, handlelength=1); ax.spines[['top', 'right']].set_visible(False); label(ax, 'a', -0.25)
ax = fig.add_axes([0.40, 0.2, 0.2, 0.68])
ax.bar(x - 0.18, SE.kw_sig, 0.34, color='#3B5B8C'); ax.bar(x + 0.18, SE.kw_study_sig, 0.34, color='#9AAFCB')
ax.set_xticks(x); ax.set_xticklabels([lab[s] for s in SE.set], fontsize=5.4); ax.set_xlabel('Genome equivalents', fontsize=6); ax.set_ylim(0, 11); ax.set_yticks(range(0, 11, 2))
for xi, (a_, b_) in enumerate(zip(SE.kw_sig, SE.kw_study_sig)): ax.text(xi - 0.18, a_ + 0.15, str(a_), ha='center', fontsize=5.4); ax.text(xi + 0.18, b_ + 0.15, str(b_), ha='center', fontsize=5.4)
ax.set_ylabel('Families differing (FDR < 0.05)'); ax.spines[['top', 'right']].set_visible(False); label(ax, 'b', -0.3)
ax = fig.add_axes([0.68, 0.2, 0.3, 0.68])
for i, f in enumerate(FAM):
    v = -np.log10(LOO[f].clip(lower=1e-300)); ax.scatter(np.full(len(v), i) + rng.uniform(-0.2, 0.2, len(v)), v, s=1.2, c='#3B5B8C', linewidths=0, rasterized=True)
ax.set_xticks(range(len(FAM))); ax.set_xticklabels(FAM, rotation=60, ha='right', rotation_mode='anchor', fontsize=5.4); ax.set_ylabel('−log$_{10}$ P, one study left out')
ax.axhline(-np.log10(0.05), color='#888888', lw=0.6, ls='--'); ax.text(len(FAM) - 0.4, -np.log10(0.05) + 0.4, 'P = 0.05', ha='right', fontsize=5.0, color='#555555')
ax.set_ylim(0, None); ax.spines[['top', 'right']].set_visible(False); label(ax, 'c', -0.2)
save(fig, f'{OUT}/Supplementary_Figure_5'); plt.close(fig)
with pd.ExcelWriter(f'{OUT}/Supplementary_Figure_5_source_data.xlsx') as xw: SE.to_excel(xw, sheet_name='ab', index=False); LOO.to_excel(xw, sheet_name='c_leave_one_study_out', index=False)

# ---- SF6 family abundance, study means
K = pd.read_csv(f'{R}/s05_family/main_kruskal.csv').set_index('family')
stm = F.groupby(['host_system', 'study'])[FAM].mean().reset_index()
fig, axs = plt.subplots(2, 5, figsize=(174 * MM, 112 * MM)); fig.subplots_adjust(left=0.08, right=0.99, bottom=0.13, top=0.95, wspace=0.38, hspace=0.62)
for ax, f in zip(axs.flat, FAM):
    for j, h in enumerate(HOSTS):
        v = stm[stm.host_system == h][f].values; strip_box(ax, np.log10(v + 1e-5), j, COL[h], 0.6, 3)
    ax.set_title(f, fontsize=6.4); ax.set_xticks(range(4)); ax.set_xticklabels(HOSTS, rotation=45, ha='right', rotation_mode='anchor', fontsize=5.2); logticks(ax, -5, 1)
    ax.text(0.5, 0.98, 'FDR < 0.001' if K.loc[f, 'FDR_study'] < 0.001 else f"FDR = {K.loc[f, 'FDR_study']:.3f}", transform=ax.transAxes, ha='center', va='top', fontsize=5.2)
    ax.set_ylim(-5.3, 1.8); ax.spines[['top', 'right']].set_visible(False)
for ax in axs[:, 0]: ax.set_ylabel('Study mean RPKG (+10$^{-5}$)')
letter_axes(axs.flat, -0.32, 1.06)
save(fig, f'{OUT}/Supplementary_Figure_6'); plt.close(fig)
stm.to_excel(f'{OUT}/Supplementary_Figure_6_source_data.xlsx', index=False)

# ---- SF7 subfamily composition per host system (all eight families with subfamily labels)
PS = pd.read_csv(f'{R}/s06_subfamily/main_subfamily_permanova.csv').set_index('family')
L = pd.read_csv(f'{R}/s04_build/subfamily_rpkg_long_main.csv'); hostof = F.set_index('sample_id').host_system
fams = list(PS.index)
fig, axs = plt.subplots(2, 4, figsize=(174 * MM, 178 * MM)); fig.subplots_adjust(left=0.07, right=0.99, bottom=0.17, top=0.95, wspace=0.3, hspace=0.95)
src = []; USED7 = []
_t = plt.get_cmap('tab20').colors; POOL7 = [_t[i] for i in [0, 2, 4, 6, 8, 10, 12, 16, 18, 1, 3, 5, 7, 9, 11, 17, 19]] + list(plt.get_cmap('tab20b').colors) + list(plt.get_cmap('Set3').colors) + list(plt.get_cmap('Dark2').colors) + list(plt.get_cmap('Paired').colors)
for ax, f in zip(axs.flat, fams):
    x = L[L.family == f].pivot_table(index='sample_id', columns='subfamily', values='RPKG', aggfunc='sum', fill_value=0)
    x = x[x.sum(axis=1) > 0]; sh = x.div(x.sum(axis=1), axis=0); sh['host'] = sh.index.map(hostof)
    m = sh.groupby('host').mean().reindex(HOSTS) * 100; nn = sh.host.value_counts().reindex(HOSTS).fillna(0).astype(int)
    cols = sorted([c for c in m.columns if c != 'unassigned'], key=lambda c: -m[c].sum())
    cols += ['unassigned'] if 'unassigned' in m.columns else []
    colors = {}
    for c in cols:
        if c == 'unassigned': colors[c] = '#F3E5AB'
        else: colors[c] = POOL7[len(USED7)]; USED7.append(c)
    bottom = np.zeros(4)
    for c in cols:
        v = m[c].fillna(0).values; ax.bar(range(4), v, bottom=bottom, color=colors[c], width=0.75, linewidth=0); bottom += v
    for k in range(4): ax.text(k, 101, str(nn.iloc[k]), ha='center', va='bottom', fontsize=4.6)
    ax.set_title(f'{f}\nR² = {PS.loc[f, "R2"]:.2f}, ' + ('P = 0.001' if PS.loc[f, 'P'] <= 0.001 else f'P = {PS.loc[f, "P"]:.3f}'), fontsize=5.8, pad=9); ax.set_xticks(range(4)); ax.set_xticklabels(HOSTS, rotation=45, ha='right', rotation_mode='anchor', fontsize=5.2)
    ax.set_ylim(0, 100); ax.spines[['top', 'right']].set_visible(False)
    ax.legend(handles=[Patch(facecolor=colors[c], label=re.sub(r'AA1\.Rdh \((\w+)\)', r'\1', c.replace('motifgroup_', 'mg '))) for c in cols],
              loc='upper left', bbox_to_anchor=(-0.15, -0.3), frameon=False, fontsize=4.4, handlelength=0.7, handletextpad=0.3, labelspacing=0.2, borderaxespad=0, ncol=2 if len(cols) > 7 else 1, columnspacing=0.8)
    src.append(m.reset_index().assign(family=f, n=nn.values))
for ax in axs[:, 0]: ax.set_ylabel('Mean share of family reads (%)')
letter_axes(axs.flat, -0.3, 1.08)
save(fig, f'{OUT}/Supplementary_Figure_7'); plt.close(fig)
pd.concat(src).to_excel(f'{OUT}/Supplementary_Figure_7_source_data.xlsx', index=False)

# ---- SF8 subfamily co-occurrence (Spearman rho, FDR < 0.05 and |rho| > 0.3)
CC = pd.read_csv(f'{R}/s06_subfamily/main_subfamily_cooccurrence.csv'); E = CC[(CC.fdr < 0.05) & (CC.rho.abs() > 0.3)].copy()
nodes = sorted(set(E.a) | set(E.b), key=lambda n: (FAM.index(n.split('::')[0]), n)); ni = {n: i for i, n in enumerate(nodes)}
Mx = np.full((len(nodes), len(nodes)), np.nan)
for r in E.itertuples(): Mx[ni[r.a], ni[r.b]] = Mx[ni[r.b], ni[r.a]] = r.rho
fig, ax = plt.subplots(figsize=(150 * MM, 140 * MM)); fig.subplots_adjust(left=0.2, right=0.92, bottom=0.2, top=0.98)
im = ax.imshow(Mx, cmap='RdBu_r', vmin=-1, vmax=1)
lab = [re.sub(r'AA1\.Rdh \((\w+)\)', r'\1', n.split('::')[1].replace('motifgroup_', 'mg ')) for n in nodes]
ax.set_xticks(range(len(nodes))); ax.set_xticklabels(lab, rotation=90, fontsize=4.6); ax.set_yticks(range(len(nodes))); ax.set_yticklabels(lab, fontsize=4.6)
prev = None
for i, n in enumerate(nodes):
    f = n.split('::')[0]
    if f != prev and i > 0: ax.axhline(i - 0.5, color='#999999', lw=0.4); ax.axvline(i - 0.5, color='#999999', lw=0.4)
    if f != prev:
        k = sum(1 for m in nodes if m.startswith(f + '::')); ax.text(-0.17, 1 - (i + k / 2) / len(nodes), f, transform=ax.transAxes, fontsize=5.2, fontweight='bold', ha='right', va='center')
    prev = f
cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02); cb.set_label('Spearman ρ', fontsize=5.6); cb.ax.tick_params(labelsize=5)
save(fig, f'{OUT}/Supplementary_Figure_8'); plt.close(fig)
E.to_excel(f'{OUT}/Supplementary_Figure_8_source_data.xlsx', index=False)

# ---- SF9 geography
GM = pd.read_csv(f'{R}/s10_geography/geo_mantel.csv')
fig, ax = plt.subplots(figsize=(84 * MM, 56 * MM)); fig.subplots_adjust(left=0.16, right=0.98, bottom=0.2, top=0.92)
for j, h in enumerate(HOSTS):
    for k, lv in enumerate(['sample', 'study']):
        r = GM[(GM.host == h) & (GM.level == lv)].iloc[0]; p = j + (k - 0.5) * 0.36
        if lv == 'sample': ax.bar(p, r.r, 0.34, color=COL[h])
        else:
            matplotlib.rcParams['hatch.linewidth'] = 0.5; ax.bar(p, r.r, 0.34, facecolor='white', edgecolor=COL[h], hatch='//////', linewidth=0.8)
        ax.text(p, max(r.r, 0) + 0.004, 'FDR ' + ('< 0.001' if r.FDR < 0.001 else f'= {r.FDR:.3f}'), ha='center', va='bottom', rotation=90, fontsize=4.8)
ax.axhline(0, color='black', lw=0.5); ax.set_xticks(range(4)); ax.set_xticklabels(HOSTS); ax.set_ylabel('Mantel r (Spearman)'); ax.set_ylim(-0.05, 0.19)
ax.legend(handles=[tup(), tup('//////')], labels=['samples', 'study means'], handler_map={tuple: HandlerTuple(ndivide=None, pad=0)}, frameon=False, fontsize=5.4, handlelength=3)
ax.spines[['top', 'right']].set_visible(False)
save(fig, f'{OUT}/Supplementary_Figure_9'); plt.close(fig)
GM.to_excel(f'{OUT}/Supplementary_Figure_9_source_data.xlsx', index=False)

# ---- SF10 host vs surroundings, both measures
AT = pd.read_csv(f'{R}/s11_ambient/host_vs_ambient_family.csv'); AT = AT[AT.set != 'NULL']
fig, axs = plt.subplots(1, 4, figsize=(174 * MM, 68 * MM), sharey=True); fig.subplots_adjust(left=0.07, right=0.99, bottom=0.3, top=0.92, wspace=0.08)
for ax, h in zip(axs, HOSTS):
    for k, (m, c_) in enumerate([('RPKG', COL[h]), ('null_normalised', '#555555')]):
        x = AT[(AT.system == h) & (AT.measure == m)].set_index('set').reindex(FAM)
        lr = np.log2(x.ratio.replace([np.inf], np.nan)); fin = np.isfinite(lr)
        pos = np.arange(len(FAM)) + (k - 0.5) * 0.38
        if m == 'RPKG': ax.bar(pos[fin], lr[fin], 0.36, color=COL[h])
        else:
            matplotlib.rcParams['hatch.linewidth'] = 0.5; ax.bar(pos[fin], lr[fin], 0.36, facecolor='white', edgecolor=COL[h], hatch='//////', linewidth=0.8)
        for i, f in enumerate(FAM):
            if np.isinf(x.ratio.iloc[i]): ax.text(pos[i], 0.3, '∞', ha='center', fontsize=5.6, color=c_)
            if x.FDR.iloc[i] < 0.05: ax.text(pos[i], (lr.iloc[i] if fin.iloc[i] else 1.2) + (0.3 if (fin.iloc[i] and lr.iloc[i] >= 0) or not fin.iloc[i] else -1.0), '*', ha='center', fontsize=6)
    ax.axhline(0, color='black', lw=0.5); ax.set_title(h, fontsize=6.4); ax.set_xticks(range(len(FAM))); ax.set_xticklabels(FAM, rotation=60, ha='right', rotation_mode='anchor', fontsize=5.2)
    ax.spines[['top', 'right']].set_visible(False)
axs[0].set_ylabel('log$_2$ host / surroundings (study means)')
letter_axes(axs, -0.12, 1.06)
fig.legend(handles=[tup(), tup('//////')], labels=['RPKG', 'null-normalised'], handler_map={tuple: HandlerTuple(ndivide=None, pad=0)}, ncol=2, frameon=False, loc='lower center', bbox_to_anchor=(0.5, 0.0), fontsize=5.6, handlelength=3)
save(fig, f'{OUT}/Supplementary_Figure_10'); plt.close(fig)
AT.to_excel(f'{OUT}/Supplementary_Figure_10_source_data.xlsx', index=False)

# ---- SF11 metatranscriptomes, study means
MS = pd.read_csv(f'{R}/s07_metatx/metatx_rpkm_per_study.csv'); MK = pd.read_csv(f'{R}/s07_metatx/metatx_family_tests.csv').set_index('family')
fig, axs = plt.subplots(2, 5, figsize=(174 * MM, 112 * MM)); fig.subplots_adjust(left=0.08, right=0.99, bottom=0.13, top=0.95, wspace=0.38, hspace=0.62)
for ax, f in zip(axs.flat, FAM):
    for j, h in enumerate(HOSTS):
        v = MS[MS.host_system == h][f].values; strip_box(ax, np.log10(v + 1e-3), j, COL[h], 0.6, 3)
    ax.set_title(f, fontsize=6.4); ax.set_xticks(range(4)); ax.set_xticklabels(HOSTS, rotation=45, ha='right', rotation_mode='anchor', fontsize=5.2); logticks(ax, -3, 3)
    q = MK.loc[f, 'FDR_study']; ax.text(0.5, 0.98, 'FDR < 0.001' if q < 0.001 else f'FDR = {q:.3f}', transform=ax.transAxes, ha='center', va='top', fontsize=5.2)
    ax.set_ylim(-3.3, 3.6); ax.spines[['top', 'right']].set_visible(False)
for ax in axs[:, 0]: ax.set_ylabel('Study mean RPKM (+10$^{-3}$)')
letter_axes(axs.flat, -0.32, 1.06)
save(fig, f'{OUT}/Supplementary_Figure_11'); plt.close(fig)
MS.to_excel(f'{OUT}/Supplementary_Figure_11_source_data.xlsx', index=False)

# ---- SF12 signal peptides
DT = pd.read_csv(f'{R}/s14_hitlevel/deeptmhmm_summary.csv'); DTT = pd.read_csv(f'{B}/inputs/hitlevel/deeptmhmm_host_tests.tsv', sep='\t').set_index('gene')
fs = [f for f in FAM if DT[DT.gene == f].signal_peptide.sum() > 0]
fig, ax = plt.subplots(figsize=(110 * MM, 60 * MM)); fig.subplots_adjust(left=0.1, right=0.99, bottom=0.16, top=0.74)
for j, h in enumerate(HOSTS):
    sub = DT[DT.environment == h].set_index('gene').reindex(fs)
    pos = np.arange(len(fs)) + (j - 1.5) * 0.2
    ax.bar(pos, sub.SP_pct.fillna(0), 0.19, color=COL[h], label=h)
    for p, n, k in zip(pos, sub.genes, sub.signal_peptide):
        if pd.notna(n): ax.text(p, 2, f'{int(k)}/{int(n)}', rotation=90, fontsize=4.4, ha='center', va='bottom')
for i_, f_ in enumerate(fs):
    q1 = DTT.loc[f_, 'FDR_raw']; q2 = DTT.loc[f_, 'FDR_order']
    ax.text(i_, 101, 'FDR ' + ('< 0.001' if q1 < 0.001 else f'{q1:.3f}') + '\nwithin order ' + ('n/a' if pd.isna(q2) else ('< 0.001' if q2 < 0.001 else f'{q2:.2f}')), ha='center', va='bottom', fontsize=4.8)
ax.set_xticks(range(len(fs))); ax.set_xticklabels(fs)
ax.set_ylabel('Genes with signal peptide (%)'); ax.set_ylim(0, 100); ax.spines[['top', 'right']].set_visible(False)
ax.legend(ncol=4, frameon=False, loc='lower center', bbox_to_anchor=(0.5, 1.2), handlelength=1)
save(fig, f'{OUT}/Supplementary_Figure_12'); plt.close(fig)
DT.to_excel(f'{OUT}/Supplementary_Figure_12_source_data.xlsx', index=False)

# ---- SF13 structural PERMANOVA + PCoA
SPm = pd.read_csv(f'{R}/s13_structures/structure_permanova_529.csv'); SPm = SPm[SPm.P.notna()].sort_values('R2')
fig = plt.figure(figsize=(174 * MM, 110 * MM))
ax = fig.add_axes([0.2, 0.56, 0.36, 0.4])
y = np.arange(len(SPm)); ax.barh(y - 0.18, SPm.R2, 0.34, color=['#0072B2' if f < 0.05 else '#F5B7A0' for f in SPm.FDR])
ax.barh(y + 0.18, SPm.R2_given_length.fillna(0), 0.34, color=['#56B4E9' if (pd.notna(f) and f < 0.05) else '#FAE0B8' for f in SPm.FDR_given_length])
ax.set_yticks(y); ax.set_yticklabels([f"{g.replace('__', ' ').replace('motifgroup_', 'mg ')} ({n})" for g, n in zip(SPm.group, SPm.n)], fontsize=5.2)
for yi, r_ in zip(y, SPm.itertuples()):
    ax.text(r_.R2 + 0.01, yi - 0.18, 'FDR ' + ('< 0.01' if r_.FDR < 0.01 else f'{r_.FDR:.2f}'), va='center', fontsize=4.4)
    if pd.notna(r_.FDR_given_length): ax.text(r_.R2_given_length + 0.01, yi + 0.18, 'FDR ' + ('< 0.01' if r_.FDR_given_length < 0.01 else f'{r_.FDR_given_length:.2f}'), va='center', fontsize=4.4)
ax.set_xlim(0, 1.0); ax.set_xlabel('R² (TM-score dissimilarity)'); ax.spines[['top', 'right']].set_visible(False)
ax.legend(handles=[Patch(facecolor='#0072B2', label='host system'), Patch(facecolor='#56B4E9', label='after chain length'), Patch(facecolor='#F5B7A0', label='FDR ≥ 0.05 (host system)'), Patch(facecolor='#FAE0B8', label='FDR ≥ 0.05 (after chain length)')], frameon=False, fontsize=5.2, loc='lower right', handlelength=1); label(ax, 'a', -0.5)
grp = ['RDase__motifgroup_CxxCxxCxxxCP', 'HyDase__motifgroup_RxxR', 'VanHase__V-BPO_ambiguous', 'VanHase__V-IPO']
for k, g in enumerate(grp):
    D = pd.read_csv(f'{R}/s13_structures/permanova_curated/{g}_dist.csv', index_col=0).values; m = pd.read_csv(f'{R}/s13_structures/permanova_curated/{g}_meta.csv')
    X, ev = pcoa(D); ax = fig.add_axes([0.08 + k * 0.245, 0.08, 0.165, 0.33])
    for h in HOSTS:
        kk = (m.host == h).values; ax.scatter(X[kk, 0], X[kk, 1], s=6, c=COL[h], linewidths=0)
    ax.set_title(g.replace('__', ' ').replace('motifgroup_', 'mg '), fontsize=5.8); ax.set_xlabel(f'PCo1 ({100 * ev[0]:.0f}%)', fontsize=5.4); ax.set_ylabel(f'PCo2 ({100 * ev[1]:.0f}%)', fontsize=5.4)
    ax.tick_params(labelsize=4.8); ax.spines[['top', 'right']].set_visible(False)
    label(ax, chr(66 + k), -0.35)
fig.legend(handles=[Patch(facecolor=COL[h], label=h) for h in HOSTS], ncol=1, frameon=False, loc='center left', bbox_to_anchor=(0.62, 0.75), handlelength=1)
save(fig, f'{OUT}/Supplementary_Figure_13'); plt.close(fig)
SPm.to_excel(f'{OUT}/Supplementary_Figure_13_source_data.xlsx', index=False)

# ---- SF14 chain length and pockets (raw, normalised and length-adjusted)
P2 = pd.read_excel(f'{B}/inputs/structures/Supplementary_Data_5.xlsx', 'p2rank_pockets')
SS = pd.read_csv(f'{R}/s13_structures/structures_529_per_structure.csv'); P2 = P2[P2.safe_id.isin(SS.safe_id)]
LA = pd.read_csv(f'{R}/s13_structures/pocket_length_adjusted_529.csv')
sets = [('RDase', 'motifgroup_CxxCxxCxxxCP'), ('HyDase', 'motifgroup_RxxR'), ('VanHase', 'V-BPO_ambiguous'), ('VanHase', 'V-IPO')]
mets = [('Chain length (residues)', 'Chain length (residues)', None),
        ('P2Rank pockets per 100 residues', 'P2Rank pockets\nper 100 residues', 'p2r100'),
        ('P2Rank mean score, top 3 pockets', 'P2Rank score,\ntop 3 pockets', 'top3')]
PT = pd.read_csv(f'{R}/s13_structures/pocket_tests_529.csv')
def ftxt(q): return 'FDR < 0.001' if q < 0.001 else f'FDR = {q:.3f}'
def raw_fdr(f, s, mtr, key):
    if key is None:
        r = PT[(PT.family == f) & (PT.subfamily == s) & (PT.metric == mtr)]
        return float(r.iloc[0].FDR) if len(r) else None
    r = LA[(LA.family == f) & (LA.subfamily == s) & (LA.metric == key)]
    return float(r.iloc[0].FDR) if len(r) else None
def adj_fdr(f, s, key):
    if key is None: return None
    r = LA[(LA.family == f) & (LA.subfamily == s) & (LA.metric == key)]
    return float(r.iloc[0].FDR_length_adjusted) if len(r) else None
fig = plt.figure(figsize=(174 * MM, 240 * MM))
gs = fig.add_gridspec(5, 12, left=0.10, right=0.99, bottom=0.02, top=0.955, wspace=2.2, hspace=0.40, height_ratios=[1, 1, 1, 1.15, 1.0])
axs = np.empty((3, 4), dtype=object)
for _i in range(3):
    for _j in range(4): axs[_i, _j] = fig.add_subplot(gs[_i, 3 * _j:3 * _j + 3])
axs_new = [fig.add_subplot(gs[3, 4 * _k:4 * _k + 4]) for _k in range(3)]
axr = fig.add_subplot(gs[4, :]); axr.imshow(plt.imread(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'assets', 'pocket_renders_4structures.png'))); axr.axis('off')
for i, (mtr, lab, key) in enumerate(mets):
    for j, (f, s) in enumerate(sets):
        ax = axs[i, j]; d = P2[(P2['Enzyme family'] == f) & (P2['Subfamily'] == s)]
        hs = [h for h in HOSTS if (d['Host system'] == h).any()]
        for k, h in enumerate(hs): strip_box(ax, d[d['Host system'] == h][mtr].values, k, COL[h], 0.6, 4)
        ax.set_xticks(range(len(hs)))
        if i == len(mets) - 1: ax.set_xticklabels(hs, rotation=45, ha='right', rotation_mode='anchor', fontsize=5.2)
        else: ax.set_xticklabels([])
        q = raw_fdr(f, s, mtr, key); qa = adj_fdr(f, s, key)
        txt = ftxt(q) if q is not None else 'not tested'
        if qa is not None: txt += '\nlength-adjusted ' + ftxt(qa)
        elif key is not None: txt += '\nlength adjustment: too few structures'
        ax.text(0.5, 1.0, txt, transform=ax.transAxes, ha='center', va='top', fontsize=4.8)
        if i == 0: ax.set_title(f'{f}\n{s.replace("motifgroup_", "mg ")}', fontsize=6)
        if j == 0: ax.set_ylabel(lab)
        ax.spines[['top', 'right']].set_visible(False); ax.margins(y=0.20)
SETC = {'GSTDase': '#E69F00', 'MetDase, NHFeHase': '#009E73'}
a = axs_new[0]
for st, c in [('RDase, HyDase, FlaHase, VanHase', '#44546A')] + list(SETC.items()):
    d = P2[P2['Structure set'] == st]
    a.scatter(d['fpocket pockets'], d['P2Rank pockets'], s=3, c=c, alpha=0.55, linewidths=0, rasterized=True, label='Main four families' if st.startswith('RDase') else st)
a.set_xlabel('Pockets per structure, fpocket'); a.set_ylabel('Pockets per structure, P2Rank')
a.text(0.04, 0.96, 'Spearman \u03c1 = 0.80\nn = 529', transform=a.transAxes, va='top', fontsize=5.6)
a.legend(frameon=False, fontsize=5, loc='lower right', handletextpad=0.2, markerscale=1.6)
b = axs_new[1]
b.scatter(P2['fpocket pockets per 100 residues'], P2['P2Rank pockets per 100 residues'], s=3, c='#44546A', alpha=0.5, linewidths=0, rasterized=True)
b.set_xlabel('fpocket pockets per 100 residues'); b.set_ylabel('P2Rank pockets per 100 residues')
b.text(0.04, 0.96, 'Spearman \u03c1 = 0.28\nn = 529', transform=b.transAxes, va='top', fontsize=5.6)
c = axs_new[2]
for fam, mk in [('RDase', 'o'), ('VanHase', 's')]:
    d = P2[(P2['Enzyme family'] == fam) & (P2['Subfamily'].isin(['motifgroup_CxxCxxCxxxCP', 'V-BPO_ambiguous']))]
    c.scatter(d['Chain length (residues)'], d['P2Rank mean score, top 3 pockets'], s=11, marker=mk, c=[COL[h] for h in d['Host system']], edgecolors='#333333', linewidths=0.3, label=fam, zorder=3)
c.set_xlabel('Chain length (residues)'); c.set_ylabel('P2Rank score, top 3 pockets')
c.legend(frameon=False, fontsize=5, loc='lower right', handletextpad=0.2, bbox_to_anchor=(1.0, 0.0))
for ax_ in axs_new: ax_.spines[['top', 'right']].set_visible(False); ax_.tick_params(length=2)
letter_axes(list(axs.flat) + axs_new, -0.3, 1.04)
axr.text(-0.012, 1.02, 'P', transform=axr.transAxes, fontsize=9, fontweight='bold', va='bottom', ha='left')
save(fig, f'{OUT}/Supplementary_Figure_14'); plt.close(fig)
with pd.ExcelWriter(f'{OUT}/Supplementary_Figure_14_source_data.xlsx') as xw:
    P2.to_excel(xw, sheet_name='pocket_metrics', index=False)
    LA.to_excel(xw, sheet_name='length_adjusted_tests', index=False)
    pd.read_csv(f'{R}/s13_structures/pocket_tool_agreement_529.csv').to_excel(xw, sheet_name='tool_agreement', index=False)
print('supplementary figures written:', sorted(os.path.basename(p) for p in glob.glob(f'{OUT}/*.pdf')))

# ---- SF15 breadth of repertoire per carrier genus
TX = pd.read_csv(f'{B}/inputs/hitlevel/All_families_taxonomy_curated9396.tsv', sep='\t')
TX['genus'] = TX.classification.fillna('').str.extract(r'g__([^;]*)')[0].fillna('')
TX = TX[TX.genus != '']
nfam = TX.groupby('genus').gene.nunique(); ngen = TX.groupby('genus').size()
order = sorted(nfam.index, key=lambda g: (-nfam[g], -ngen[g]))[:20]
fig, axs = plt.subplots(1, 2, figsize=(174 * MM, 92 * MM), gridspec_kw=dict(width_ratios=[1, 1.45], wspace=0.06))
ax = axs[0]; left = np.zeros(len(order))
for h in HOSTS:
    v = np.array([((TX.genus == g) & (TX.environment == h)).sum() for g in order], float)
    ax.barh(range(len(order)), v, left=left, color=COL[h], height=0.72, lw=0.3, edgecolor='white', label=h); left += v
ax.set_yticks(range(len(order)))
ax.set_yticklabels(order, fontsize=5.4, fontstyle='italic')
for lbl in ax.get_yticklabels():
    if any(ch.isdigit() for ch in lbl.get_text()) or lbl.get_text().isupper(): lbl.set_fontstyle('normal')
ax.invert_yaxis(); ax.set_xlabel('Curated genes (n)'); ax.set_xscale('log'); ax.set_xlim(0.8, max(left) * 1.6)
ax.spines[['top', 'right']].set_visible(False); ax.tick_params(length=2)
ax.legend(frameon=False, fontsize=5.2, loc='lower right', handlelength=1.0, handleheight=0.8, borderpad=0.2)
bx = axs[1]
for i, g in enumerate(order):
    sub = TX[TX.genus == g]
    for j, f in enumerate(FAM):
        d = sub[sub.gene == f]
        if len(d):
            h = d.environment.value_counts().idxmax()
            bx.scatter(j, i, s=5 + 9 * np.sqrt(len(d)), c=[COL[h]], edgecolors='#555555', linewidths=0.3, zorder=3)
    bx.text(len(FAM) - 0.2, i, str(int(nfam[g])), va='center', ha='left', fontsize=5.4, fontweight='bold')
bx.set_xlim(-0.6, len(FAM) + 0.5); bx.set_ylim(len(order) - 0.5, -0.5)
bx.set_xticks(range(len(FAM))); bx.set_xticklabels(FAM, rotation=90, fontsize=5.4)
bx.set_yticks(range(len(order))); bx.set_yticklabels([])
bx.text(len(FAM) - 0.2, -0.9, 'families', fontsize=5.4, fontweight='bold', ha='left', va='bottom')
bx.grid(color='#EEEEEE', lw=0.4); bx.set_axisbelow(True); bx.tick_params(length=0)
for s in ['top', 'right', 'left', 'bottom']: bx.spines[s].set_visible(False)
for k_, n_ in enumerate([1, 10, 100]):
    kx = fig.add_axes([0.70 + k_ * 0.055, 0.028, 0.02, 0.045]); kx.axis('off'); kx.set_xlim(0, 1); kx.set_ylim(0, 1)
    kx.scatter([0.5], [0.5], s=5 + 9 * np.sqrt(n_), c='#BBBBBB', edgecolors='#555555', linewidths=0.3, clip_on=False)
    kx.text(1.6, 0.5, str(n_), va='center', fontsize=5.2)
fig.text(0.625, 0.050, 'Curated genes', fontsize=5.6, fontweight='bold', va='center', ha='right')
letter_axes(axs, -0.22, 1.02)
fig.subplots_adjust(left=0.135, right=0.975, top=0.93, bottom=0.20)
save(fig, f'{OUT}/Supplementary_Figure_15'); plt.close(fig)
with pd.ExcelWriter(f'{OUT}/Supplementary_Figure_15_source_data.xlsx') as xw:
    pd.DataFrame({'genus': order, 'families': [int(nfam[g]) for g in order], 'curated_genes': [int(ngen[g]) for g in order]}).to_excel(xw, sheet_name='top20_genera', index=False)
    TX.groupby(['genus', 'gene', 'environment']).size().rename('genes').reset_index().to_excel(xw, sheet_name='genus_family_host', index=False)
print('SF15 written; Streptomyces families =', int(nfam['Streptomyces']))

# ---- publication numbering: figures are numbered by first citation in the manuscript (old -> new); done once, after generation
import glob as _g
FMAP = {7: 10, 8: 11, 9: 7, 10: 8, 11: 9}
_tmp = []
for _o, _n in FMAP.items():
    for _f in _g.glob(f'{OUT}/Supplementary_Figure_{_o}[._]*'):
        _t = _f.replace(f'Supplementary_Figure_{_o}', f'TMPFIG_{_n}'); os.rename(_f, _t); _tmp.append(_t)
for _t in _tmp: os.rename(_t, _t.replace('TMPFIG_', 'Supplementary_Figure_'))
print('Supplementary figures renumbered:', FMAP)
