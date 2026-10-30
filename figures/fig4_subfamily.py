"""Figure 4. Subfamilies partition among host systems more strongly than families.
usage: python3 fig4_subfamily.py <FINAL dir>
a variance in composition explained by host system: ten-family composition vs subfamily composition within each family
  (sample and study level; PERMANOVA, Bray-Curtis); b all 48 indicator subfamilies (multipatt IndVal, labels in grey rest on <= 3 curated genes;
  FDR <= 0.05): dots = share of samples with reads, filled = host system(s) of the indicator combination, with IndVal, curated
  genes and the genome-level lineage control (GTDB order); c mean subfamily composition per host system (samples with reads
  on the family); e three lines of host support for the six subfamilies that stay host-associated within GTDB order:
  within-order odds ratio, null-normalized fold above surrounding metagenomes (s16), and share of metatranscriptome runs with reads.
"""
import sys, os, re, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import *
from matplotlib.patches import Patch
B = sys.argv[1]; T = f'{B}/figures/tmp_fig4'; os.makedirs(T, exist_ok=True); out = f'{B}/figures/Figure_4'
V = pd.read_csv(f'{B}/results/s05_family/main_vegan_tests.csv').set_index('test')
PS = pd.read_csv(f'{B}/results/s06_subfamily/main_subfamily_permanova.csv').set_index('family')
I = pd.read_csv(f'{B}/results/s06_subfamily/main_subfamily_indicators_significant.csv')
L = pd.read_csv(f'{B}/results/s04_build/subfamily_rpkg_long_main.csv')
F = pd.read_csv(f'{B}/results/s04_build/family_rpkg_main.csv')
host = F.set_index('sample_id').host_system
KS = pd.read_csv(f'{B}/results/s06_subfamily/main_subfamily_kruskal.csv')


# ---- p: combined PCoA of all labelled subfamilies
from scipy.stats import chi2
from matplotlib.patches import Ellipse
PC = pd.read_csv(f'{B}/results/s06_subfamily/combined_pcoa_coordinates.csv'); CP = pd.read_csv(f'{B}/results/s06_subfamily/combined_permanova.csv').set_index('level')
fig, ax = plt.subplots(figsize=(84 * MM, 62 * MM)); fig.subplots_adjust(left=0.13, right=0.98, bottom=0.17, top=0.82)
rng_ = np.random.default_rng(2)
for h in HOSTS:
    x = PC[PC.host_system == h]
    ax.scatter(x.PCo1, x.PCo2, s=5, c=COL[h], alpha=0.75, edgecolors='white', linewidths=0.15, label=f'{h} (n = {len(x)})')
for h in HOSTS:
    x = PC[PC.host_system == h]; c = np.cov(x.PCo1, x.PCo2); w, v_ = np.linalg.eigh(c); k = np.sqrt(chi2.ppf(0.95, 2))
    ax.add_patch(Ellipse((x.PCo1.mean(), x.PCo2.mean()), 2 * k * np.sqrt(w[1]), 2 * k * np.sqrt(w[0]), angle=np.degrees(np.arctan2(v_[1, 1], v_[0, 1])), fill=False, ec=COL[h], lw=1.0))
ax.set_xlabel(f'PCo1 ({100 * CP.loc["sample", "PCo1_var"]:.1f}%)'); ax.set_ylabel(f'PCo2 ({100 * CP.loc["sample", "PCo2_var"]:.1f}%)'); ax.spines[['top', 'right']].set_visible(False)
ax.set_ylim(-0.48, 1.0)
ax.legend(loc='lower center', bbox_to_anchor=(0.5, 1.0), ncol=2, frameon=False, markerscale=1.8, handletextpad=0.2, columnspacing=1.0, borderaxespad=0.2)
ax.text(0.98, 0.99, f"PERMANOVA R² = {CP.loc['sample', 'R2']:.3f}, P = {CP.loc['sample', 'P']:.3f}\nPERMDISP F = {CP.loc['sample', 'permdisp_F']:.1f}, P = {CP.loc['sample', 'permdisp_P']:.3f}", transform=ax.transAxes, ha='right', va='top', fontsize=5.4)
fig.savefig(f'{T}/p.pdf'); plt.close(fig)

# ---- a: R2
fams = [f for f in FAM if f in PS.index]
fig, ax = plt.subplots(figsize=(84 * MM, 62 * MM)); fig.subplots_adjust(left=0.2, right=0.97, bottom=0.17, top=0.88)
y = np.arange(len(fams) + 1)
r2s = [V.loc['host_sample', 'R2']] + [PS.loc[f, 'R2'] for f in fams]; r2t = [V.loc['host_study', 'R2']] + [PS.loc[f, 'R2_study'] for f in fams]
pst = [V.loc['host_study', 'P']] + [PS.loc[f, 'P_study'] for f in fams]
ax.barh(y - 0.19, r2s, height=0.36, color='#3B5B8C', label='sample level')
ax.barh(y + 0.19, r2t, height=0.36, color='#9AAFCB', label='study level')
for yi, v, p in zip(y, r2t, pst):
    if p > 0.05: ax.text(v + 0.006, yi + 0.19, 'n.s.', va='center', fontsize=5)
ax.set_yticks(y); ax.set_yticklabels(['All families'] + [f'{f} ({int(PS.loc[f, "n_subfamilies"])})' for f in fams])
ax.get_yticklabels()[0].set_fontweight('bold'); ax.invert_yaxis(); ax.axvline(V.loc['host_sample', 'R2'], color='#3B5B8C', lw=0.5, ls=':')
ax.set_xlabel('Variance explained by host system (R²)'); ax.set_xlim(0, 0.55); ax.spines[['top', 'right']].set_visible(False)
ax.legend(frameon=False, loc='lower center', bbox_to_anchor=(0.45, 1.0), ncol=2, handlelength=1)
fig.savefig(f'{T}/a.pdf'); plt.close(fig)

# ---- b: indicator subfamilies with >= 4 curated genes
J = I.copy(); J['fo'] = J.family.map({f: i for i, f in enumerate(FAM)}); J = J.sort_values(['fo', 'stat'], ascending=[True, False])
LV = {'host-associated within lineage': ('host-linked', '#1B7837'), 'not distinguishable from lineage': ('lineage-linked', '#B35806'),
      'not testable within lineage': ('not testable', '#777777'), 'too few carrier genomes': ('too few genomes', '#777777')}
n = len(J); rowh = 2.6
fig = plt.figure(figsize=(88 * MM, (n * rowh + 26) * MM)); ax = fig.add_axes([0.38, 12 / (n * rowh + 26), 0.6, n * rowh / (n * rowh + 26)])
yy = np.arange(n)
for j, h in enumerate(HOSTS):
    pv = J[f'prev_{h}'].values * 100; inc = J.combination.str.split('+').apply(lambda c: h in c).values
    ax.scatter(np.full(n, j)[inc], yy[inc], s=3 + pv[inc] * 0.75, c=COL[h], edgecolors='none', zorder=3)
    ax.scatter(np.full(n, j)[~inc], yy[~inc], s=3 + pv[~inc] * 0.75, facecolors='none', edgecolors='#AAAAAA', linewidths=0.5, zorder=3)
ax.set_xlim(-0.6, 8.8); ax.set_ylim(n - 0.5, -0.5); ax.set_xticks(range(4)); ax.set_xticklabels(HOSTS, rotation=40, ha='left', rotation_mode='anchor'); ax.xaxis.tick_top()
ax.set_yticks(yy); ax.set_yticklabels([re.sub(r'AA1\.Rdh \((\w+)\)', r'\1', s.replace('motifgroup_', 'mg ')) for s in J.subfamily], fontsize=5.0)
for t_, g_ in zip(ax.get_yticklabels(), J.n_genes): t_.set_color('#222222' if g_ >= 4 else '#999999')
ax.tick_params(axis='both', length=0); [s_.set_visible(False) for s_ in ax.spines.values()]
prev_f = None
for i, r in enumerate(J.itertuples()):
    if r.family != prev_f:
        if i > 0: ax.axhline(i - 0.5, color='#CCCCCC', lw=0.5, xmin=-0.95, xmax=1.0, clip_on=False)
        k = (J.family == r.family).sum(); ax.text(-0.33 - 0.0, i + (k - 1) / 2, '', fontsize=5)
        fig.text(0.01, ax.transData.transform((0, i + (k - 1) / 2))[1] / fig.dpi / fig.get_size_inches()[1], r.family, fontsize=5.4, fontweight='bold', va='center')
        prev_f = r.family
    ax.text(3.9, i, f'{r.stat:.2f}', va='center', fontsize=5.0, color='#222222' if r.n_genes >= 4 else '#999999')
    ax.text(6.0, i, f'{r.n_genes}', va='center', fontsize=5.0, ha='right', color='#222222' if r.n_genes >= 4 else '#999999')
    lab, c = LV[r.lineage_order_verdict]; ax.text(6.3, i, lab, va='center', fontsize=4.8, color=c)
for xx, t, ha in [(3.9, 'IndVal', 'left'), (6.0, 'genes', 'right'), (6.3, 'Within-order\ntest', 'left')]:
    ax.text(xx, -1.2, t, fontsize=5.4, fontweight='bold', ha=ha, va='bottom')
ax.text(-0.5, n + 1.6, 'Samples with reads (%)', fontsize=5.2, fontweight='bold', va='center', ha='left', clip_on=False)
for k_, pv in enumerate([10, 50, 100]):
    ax.scatter(-0.2 + k_ * 1.6, n + 3.4, s=3 + pv * 0.75, c='#999999', clip_on=False); ax.text(0.25 + k_ * 1.6, n + 3.4, f'{pv}', va='center', fontsize=5.2)
ax.scatter(5.0, n + 1.6, s=14, c='#666666', clip_on=False); ax.text(5.3, n + 1.6, 'indicator host', va='center', fontsize=5.2)
ax.scatter(5.0, n + 3.4, s=14, facecolors='none', edgecolors='#AAAAAA', linewidths=0.5, clip_on=False); ax.text(5.3, n + 3.4, 'other host', va='center', fontsize=5.2)

fig.savefig(f'{T}/b.pdf'); plt.close(fig)
hb = n * rowh + 26

# ---- c: subfamily composition (mean share among samples with reads on the family)
show = ['RDase', 'GSTDase', 'MetDase', 'VanHase']
fig, axs = plt.subplots(1, 4, figsize=(84 * MM, 72 * MM)); fig.subplots_adjust(left=0.13, right=0.99, bottom=0.45, top=0.91, wspace=0.3)
comp_rows = []; USED = []
_t = plt.get_cmap('tab20').colors; POOL = [_t[i] for i in [0, 2, 4, 6, 8, 10, 12, 16, 18, 1, 3, 5, 7, 9, 11, 17, 19]] + list(plt.get_cmap('tab20b').colors)
for ax, f in zip(axs, show):
    x = L[L.family == f].pivot_table(index='sample_id', columns='subfamily', values='RPKG', aggfunc='sum', fill_value=0)
    x = x[x.sum(axis=1) > 0]; sh = x.div(x.sum(axis=1), axis=0); sh['host'] = sh.index.map(host)
    m = sh.groupby('host').mean().reindex(HOSTS) * 100
    cols = [c for c in m.columns if c != 'unassigned']; cols = sorted(cols, key=lambda c: -m[c].sum())
    cols += ['unassigned'] if 'unassigned' in m.columns else []
    colors = {}
    for c in cols:
        if c == 'unassigned': colors[c] = '#F3E5AB'
        else: colors[c] = POOL[len(USED)]; USED.append(c)
    bottom = np.zeros(4)
    for c in cols:
        v = m[c].fillna(0).values; ax.bar(range(4), v, bottom=bottom, color=colors[c], width=0.75, linewidth=0, label=c.replace('motifgroup_', 'mg ')); bottom += v
    nn = sh.host.value_counts().reindex(HOSTS).fillna(0).astype(int)
    for j in range(4):
        ax.text(j, 101, f'{nn.iloc[j]}', fontsize=4.4, ha='center', va='bottom')
    ax.set_title(f, fontsize=6.2, pad=7); ax.set_xticks(range(4)); ax.set_xticklabels(HOSTS, rotation=90, fontsize=5.2)
    ax.set_ylim(0, 100); ax.spines[['top', 'right']].set_visible(False)
    if f == 'RDase': ax.set_ylabel('Mean share of family reads (%)'); ax.text(-0.7, 101, 'n', fontsize=4.4, ha='right', va='bottom')
    else: ax.set_yticklabels([])
    leg = [Patch(facecolor=colors[c], label=re.sub(r'AA1\.Rdh \((\w+)\)', r'\1', c.replace('motifgroup_', 'mg '))) for c in cols]
    ax.legend(handles=leg, loc='upper left', bbox_to_anchor=(-0.2, -0.5), frameon=False, fontsize=4.4, handlelength=0.7, handletextpad=0.3, labelspacing=0.2, borderaxespad=0)
    m.reset_index().assign(family=f).pipe(comp_rows.append)
fig.savefig(f'{T}/c.pdf'); plt.close(fig)

# ---- e: three lines of host support for the subfamilies that stay host-associated within GTDB order (s16)
SV = pd.read_csv(f'{B}/results/s16_indicator_support/indicator_host_vs_surroundings.csv')
EX = pd.read_csv(f'{B}/results/s16_indicator_support/indicator_expression.csv')
IA = pd.read_csv(f'{B}/results/s06_subfamily/main_subfamily_indicators_significant.csv')
six = [('GSTDase', 'motifgroup_SxC'), ('GSTDase', 'linD'), ('HyDase', 'dhmA'), ('HyDase', 'motifgroup_DRG'), ('VanHase', 'V-HPO'), ('VanHase', 'V-IPO_ambiguous')]
fig, axs = plt.subplots(1, 3, figsize=(84 * MM, 58 * MM), sharey=True)
fig.subplots_adjust(left=0.27, right=0.985, bottom=0.36, top=0.83, wspace=0.16)
ny = len(six); drows = []
CAPX = [3000, 300]
for i, (f, sf) in enumerate(six):
    y = ny - 1 - i
    ia = IA[(IA.family == f) & (IA.subfamily == sf)].iloc[0]; hosts = ia.combination.split('+')
    hc = COL[hosts[0]] if len(hosts) == 1 else '#555555'
    # (1) within-order odds ratio
    orr = ia.lineage_order_OR; xo = CAPX[0] if not np.isfinite(orr) else orr
    axs[0].scatter([xo], [y], s=14, marker='>' if not np.isfinite(orr) else 'o', color=hc, zorder=3, linewidths=0)
    axs[0].text(xo * (0.55 if not np.isfinite(orr) else 1.5), y + 0.28, f'{int(ia.lineage_order_strata)}', fontsize=4.6, ha='center', va='bottom', color='#333333')
    # (2) above surroundings, null-normalised fold, per indicator host
    for j, h in enumerate(hosts):
        r = SV[(SV.family == f) & (SV.subfamily == sf) & (SV.system == h)].iloc[0]
        v = r.ratio_norm; xv = CAPX[1] if not np.isfinite(v) else v; yy = y + (0.17 if j == 0 else -0.17) * (len(hosts) > 1)
        filled = r.call == 'host-enriched'
        axs[1].scatter([xv], [yy], s=14, marker='>' if not np.isfinite(v) else 'o', facecolors=COL[h] if filled else 'white', edgecolors=COL[h], linewidths=0.7, zorder=3)
        drows.append(dict(family=f, subfamily=sf, system=h, evidence='surroundings', value=v, call=r.call, FDR_RPKG=r.FDR_RPKG, FDR_null_normalised=r.FDR_norm,
                          host_studies=r.n_host_studies_RPKG, surrounding_studies=r.n_ambient_studies_RPKG))
    # (3) expression: % of runs with reads, indicator host(s) vs the other host systems pooled
    e = EX[(EX.family == f) & (EX.subfamily == sf)].iloc[0]
    oth = [h for h in HOSTS if h not in hosts]
    po = 100 * sum(e[f'det_runs_{h}'] * e[f'runs_{h}'] for h in oth) / sum(e[f'runs_{h}'] for h in oth)
    axs[2].scatter([po], [y], s=10, facecolors='white', edgecolors='#9a9a9a', linewidths=0.6, zorder=2)
    filled = (e.FDR_study < 0.05) and bool(e.top_in_indicator)
    for j, h in enumerate(hosts):
        yy = y + (0.17 if j == 0 else -0.17) * (len(hosts) > 1)
        axs[2].scatter([100 * e[f'det_runs_{h}']], [yy], s=14, facecolors=COL[h] if filled else 'white', edgecolors=COL[h], linewidths=0.7, zorder=3)
    drows.append(dict(family=f, subfamily=sf, system='+'.join(hosts), evidence='lineage (GTDB order)', value=orr, strata=ia.lineage_order_strata, FDR=ia.lineage_order_FDR))
    drows.append(dict(family=f, subfamily=sf, system='+'.join(hosts), evidence='expression', other_hosts_pct_runs=po, FDR_run=e.FDR_run, FDR_study=e.FDR_study,
                      **{f'pct_runs_{h}': 100 * e[f'det_runs_{h}'] for h in HOSTS}))
lab = [f'{f} {sf.replace("motifgroup_", "mg ").replace("_ambiguous", "_amb.")}' for f, sf in six]
axs[0].set_yticks(range(ny)); axs[0].set_yticklabels(lab[::-1], fontsize=5.4); axs[0].set_ylim(-0.6, ny - 0.4)
for ax in axs:
    ax.spines[['top', 'right']].set_visible(False); ax.tick_params(axis='x', labelsize=5.2)
    for k in range(ny): ax.axhline(k, color='#eeeeee', lw=0.4, zorder=0)
axs[0].set_xscale('log'); axs[0].set_xlim(1, 6000); axs[0].axvline(1, color='#888888', lw=0.5, ls=':')
axs[0].set_xticks([1, 10, 100, 1000]); axs[0].set_xticklabels(['1', '10', '100', '10$^3$'])
axs[0].set_xlabel('Odds ratio', fontsize=5.6); axs[0].set_title('Within bacterial\norder', fontsize=5.8, pad=4)
axs[1].set_xscale('log'); axs[1].set_xlim(0.3, 600); axs[1].axvline(1, color='#888888', lw=0.5, ls=':')
axs[1].set_xticks([1, 10, 100]); axs[1].set_xticklabels(['1', '10', '100'])
axs[1].set_xlabel('Fold, null-normalized', fontsize=5.6); axs[1].set_title('Above\nsurroundings', fontsize=5.8, pad=4)
axs[2].set_xlim(-3, 60); axs[2].set_xlabel('Runs with reads (%)', fontsize=5.6); axs[2].set_title('Transcribed in\nmetatranscriptomes', fontsize=5.8, pad=4)
from matplotlib.lines import Line2D
hl = [Line2D([], [], ls='', marker='o', ms=3.4, mfc=COL[h], mec=COL[h], label=h) for h in HOSTS]
hl += [Line2D([], [], ls='', marker='o', ms=3.4, mfc='#555555', mec='#555555', label='Coral + Rhizosphere'),
       Line2D([], [], ls='', marker='o', ms=3.4, mfc='white', mec='#555555', mew=0.7, label='FDR > 0.05'),
       Line2D([], [], ls='', marker='>', ms=3.4, mfc='#555555', mec='#555555', label='zero in comparison'),
       Line2D([], [], ls='', marker='o', ms=3.0, mfc='white', mec='#9a9a9a', mew=0.6, label='other host systems')]
for ax in axs: ax.minorticks_off()
fig.legend(handles=hl, loc='lower center', ncol=4, frameon=False, fontsize=4.6, handletextpad=0.2, columnspacing=0.7, bbox_to_anchor=(0.52, 0.0))
fig.savefig(f'{T}/d.pdf'); plt.close(fig)

H = assemble([('a', f'{T}/p.pdf', 0, 0, 84, 62), ('b', f'{T}/a.pdf', 90, 0, 84, 62), ('c', f'{T}/b.pdf', 0, 65, 88, hb),
              ('d', f'{T}/c.pdf', 90, 65, 84, 72), ('e', f'{T}/d.pdf', 90, 139, 84, 58)], out)
print('Figure 4 height mm', round(H, 1), '| indicators shown', n, 'of', len(I))
assert H <= 225, H
with pd.ExcelWriter(out + '_source_data.xlsx') as xw:
    PS.reset_index().to_excel(xw, sheet_name='b_subfamily_permanova', index=False); V.reset_index().to_excel(xw, sheet_name='b_family_permanova', index=False)
    I.to_excel(xw, sheet_name='c_all_48_indicators', index=False)
    pd.concat(comp_rows).to_excel(xw, sheet_name='d_mean_share_pct', index=False)
    pd.DataFrame(drows).to_excel(xw, sheet_name='e_host_support_six', index=False); PC.to_excel(xw, sheet_name='a_combined_pcoa', index=False); CP.reset_index().to_excel(xw, sheet_name='a_combined_permanova', index=False); KS.to_excel(xw, sheet_name='e_kruskal_all_subfamilies', index=False)
