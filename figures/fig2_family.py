"""Figure 2. Host system structures the ten-family repertoire beyond study effects (main set, 1,787 samples).
usage: python3 fig2_family.py <FINAL dir>
Panels: a NMDS (Bray-Curtis on RPKG of the ten families; samples with >= 1 target read) with 95% covariance ellipses;
b variance explained by host system (sample level, study level, study nested in host, host after log10 genome equivalents);
c indicator families (IndVal.g); d Shannon diversity of family composition with Dunn letters; e mean relative abundance (%);
f mean RPKG per family and host (log scale) with Dunn letters from the sample-level Kruskal-Wallis tests.
Source data for every panel are written next to the figure (Fig2_source_data.xlsx).
"""
import sys, os, json, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from figstyle import *
from stats_utils import kw, dunn, cld, shannon
from scipy.stats import chi2
from matplotlib.patches import Ellipse
B = sys.argv[1]; R = f'{B}/results/s05_family'; out = f'{B}/figures/Figure_2'
d = pd.read_csv(f'{B}/results/s04_build/family_rpkg_main.csv')
nm = pd.read_csv(f'{R}/main_vegan_nmds.csv'); V = pd.read_csv(f'{R}/main_vegan_tests.csv').set_index('test')
G = pd.read_csv(f'{R}/main_vegan_logGE_first.csv').set_index('term'); I = pd.read_csv(f'{R}/main_indicator_families.csv')
K = pd.read_csv(f'{R}/main_kruskal.csv').set_index('family'); RA = pd.read_csv(f'{R}/main_relative_abundance.csv')
S = json.load(open(f'{R}/main_summary.json'))
fig = plt.figure(figsize=(174 * MM, 222 * MM))
# a NMDS
ax = fig.add_axes([0.08, 0.675, 0.41, 0.265])
for h in HOSTS:
    x = nm[nm.group == h]
    ax.scatter(x.NMDS1, x.NMDS2, s=5.5, c=COL[h], alpha=0.75, edgecolors='white', linewidths=0.15, label=f'{h} (n = {len(x)})')
    c = np.cov(x.NMDS1, x.NMDS2); w, v = np.linalg.eigh(c); k = np.sqrt(chi2.ppf(0.95, 2))
    ang = np.degrees(np.arctan2(v[1, 1], v[0, 1]))
    ax.add_patch(Ellipse((x.NMDS1.mean(), x.NMDS2.mean()), 2 * k * np.sqrt(w[1]), 2 * k * np.sqrt(w[0]), angle=ang, fill=False, ec=COL[h], lw=1.0))
ax.set_xlabel('NMDS1'); ax.set_ylabel('NMDS2'); ax.spines[['top', 'right']].set_visible(False)
ax.legend(loc='lower center', bbox_to_anchor=(0.5, 1.0), ncol=2, frameon=False, markerscale=1.8, handletextpad=0.2, columnspacing=1.0, borderaxespad=0.2)
ax.text(0.02, 0.98, f"Stress = {nm.stress.iat[0]:.3f}\nPERMANOVA R² = {V.loc['host_sample', 'R2']:.3f}, P = {V.loc['host_sample', 'P']:.3f}\n"
        f"PERMDISP F = {V.loc['permdisp_sample', 'F']:.1f}, P = {V.loc['permdisp_sample', 'P']:.3f}", transform=ax.transAxes, va='top', fontsize=5.8)
label(ax, 'a', -0.17, 1.11)
# b R2
ax = fig.add_axes([0.61, 0.675, 0.37, 0.265])
bars = [('Host system\nsample level', V.loc['host_sample', 'R2'], V.loc['host_sample', 'P'], f"n = {int(V.loc['host_sample', 'n']):,} samples", '#1F4E79'),
        ('Host system\nstudy level', V.loc['host_study', 'R2'], V.loc['host_study', 'P'], f"n = {int(V.loc['host_study', 'n'])} studies", '#2E86C1'),
        ('Study nested\nin host system', V.loc['block_nested', 'R2'], V.loc['block_nested', 'P'], f"{int(V.loc['block_nested', 'df']) + 4} studies", '#7D3C98'),
        ('Host system after\ngenome equivalents', G.loc['host_system', 'R2'], G.loc['host_system', 'P'], f"GE R² = {G.loc['log10(GE)', 'R2']:.3f}", '#17A589')]
for i, (lab, r2, p, n, c) in enumerate(bars):
    ax.bar(i, r2, color=c, width=0.65); ax.text(i, r2 + 0.006, f'R² = {r2:.3f}\nP = {p:.4f}'.replace('P = 0.0010', 'P = 0.001') + f'\n{n}', ha='center', va='bottom', fontsize=5.4)
ax.set_xticks(range(4)); ax.set_xticklabels([b[0] for b in bars], fontsize=5.6); ax.set_ylim(0, 0.45); ax.set_ylabel('Variance explained (R²)')
ax.spines[['top', 'right']].set_visible(False); label(ax, 'b', -0.2, 1.11)
# c indicator families
ax = fig.add_axes([0.10, 0.37, 0.30, 0.25])
M = np.full((len(FAM), 4), np.nan)
for i, f in enumerate(FAM):
    r = I.set_index('family').loc[f]
    if pd.notna(r.FDR) and r.FDR <= 0.05:
        for h in r.combination.split('+'): M[i, HOSTS.index(h)] = r.stat
im = ax.imshow(M, cmap='Blues', vmin=0.4, vmax=1.0, aspect='auto')
for i in range(len(FAM)):
    for j in range(4):
        if not np.isnan(M[i, j]): ax.text(j, i, f'{M[i, j]:.2f}', ha='center', va='center', fontsize=5.6, color='white' if M[i, j] > 0.75 else 'black')
ax.set_xticks(range(4)); ax.set_xticklabels(HOSTS, rotation=30, ha='right'); ax.set_yticks(range(len(FAM))); ax.set_yticklabels(FAM)
ax.set_xticks(np.arange(-.5, 4), minor=True); ax.set_yticks(np.arange(-.5, len(FAM)), minor=True); ax.grid(which='minor', color='#DDDDDD', lw=0.4); ax.tick_params(which='minor', length=0)
cb = fig.colorbar(im, ax=ax, fraction=0.06, pad=0.03); cb.set_label('IndVal', fontsize=6); cb.ax.tick_params(labelsize=5.5)
label(ax, 'c', -0.30)
# d Shannon
ax = fig.add_axes([0.60, 0.37, 0.38, 0.25])
dn = d[d[FAM].sum(axis=1) > 0].copy(); dn['sh'] = shannon(dn[FAM].values); g = dn.host_system.values
H, p = kw(dn.sh.values, g); pr = dunn(dn.sh.values, g); rk = pd.Series(dn.sh).rank().values
let = cld(HOSTS, pr, {h: rk[g == h].mean() for h in HOSTS})
bp = ax.boxplot([dn.loc[g == h, 'sh'] for h in HOSTS], widths=0.6, patch_artist=True, showfliers=False, medianprops=dict(color='black', lw=0.8), whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6))
for b_, h in zip(bp['boxes'], HOSTS): b_.set_facecolor(COL[h]); b_.set_alpha(0.85); b_.set_linewidth(0.6)
rng = np.random.default_rng(1)
for i, h in enumerate(HOSTS):
    y = dn.loc[g == h, 'sh'].values; dark = tuple(0.55 * c_ for c_ in matplotlib.colors.to_rgb(COL[h]))
    ax.scatter(i + 1 + rng.uniform(-0.2, 0.2, len(y)), y, s=2.2, c=[dark], alpha=0.5, linewidths=0)
    ax.text(i + 1, 2.38, let[h], ha='center', fontweight='bold', fontsize=6.5)
ax.set_ylim(-0.05, 2.5); ax.set_xticks(range(1, 5)); ax.set_xticklabels([f'{h}\n(n = {(g == h).sum()})' for h in HOSTS]); ax.set_ylabel('Shannon diversity (families)')
ax.spines[['top', 'right']].set_visible(False); ax.text(0.5, -0.20, f'Kruskal–Wallis H = {H:.1f}, {pfmt(p)}', transform=ax.transAxes, ha='center', fontsize=5.6)
label(ax, 'd', -0.15)
# e relative abundance
ax = fig.add_axes([0.10, 0.05, 0.30, 0.25])
RAp = RA.pivot(index='family', columns='host', values='mean').reindex(index=FAM, columns=HOSTS)
order = RAp.mean(axis=1).sort_values(ascending=False).index; RAp = RAp.loc[order]
im = ax.imshow(RAp.values, cmap='Reds', aspect='auto', norm=matplotlib.colors.PowerNorm(0.5, vmin=0, vmax=max(70, RAp.values.max())))
for i in range(RAp.shape[0]):
    for j in range(4):
        v = RAp.values[i, j]; ax.text(j, i, f'{v:.1f}' if v < 10 else f'{v:.0f}', ha='center', va='center', fontsize=5.6, color='white' if v > 30 else 'black')
nre = RA.groupby('host').n.first()
ax.set_xticks(range(4)); ax.set_xticklabels([f'{h}\n(n = {nre[h]})' for h in HOSTS], fontsize=5.6); ax.set_yticks(range(len(order))); ax.set_yticklabels(order)
cb = fig.colorbar(im, ax=ax, fraction=0.06, pad=0.03); cb.set_label('Mean relative abundance (%)', fontsize=6); cb.ax.tick_params(labelsize=5.5)
label(ax, 'e', -0.30)
# f mean RPKG with letters
ax = fig.add_axes([0.60, 0.07, 0.38, 0.23])
w = 0.2
for j, h in enumerate(HOSTS):
    m = [K.loc[f, f'mean_{h}'] for f in FAM]
    ax.bar(np.arange(len(FAM)) + (j - 1.5) * w, np.array(m) + 1e-6, width=w, color=COL[h], label=h)
    for i, f in enumerate(FAM): ax.text(i + (j - 1.5) * w, max(m[i], 1e-5) * 1.5, K.loc[f, f'letter_{h}'], ha='left', va='center', fontsize=5.0, rotation=90, rotation_mode='anchor', fontweight='bold')
ax.set_yscale('log'); ax.set_ylim(1e-5, 600); ax.set_xticks(range(len(FAM))); ax.set_xticklabels(FAM, rotation=45, ha='right', rotation_mode='anchor')
ax.set_ylabel('Mean RPKG'); ax.yaxis.set_major_locator(matplotlib.ticker.LogLocator(base=10, numticks=12)); ax.yaxis.set_minor_locator(matplotlib.ticker.LogLocator(base=10, subs=np.arange(2, 10), numticks=12)); ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter()); ax.spines[['top', 'right']].set_visible(False)
ax.legend(ncol=2, frameon=False, loc='upper right', handlelength=1, columnspacing=0.8, fontsize=5.6)
label(ax, 'f', -0.15)
os.makedirs(f'{B}/figures', exist_ok=True); save(fig, out)
with pd.ExcelWriter(out + '_source_data.xlsx') as xw:
    nm.to_excel(xw, sheet_name='a_NMDS', index=False); pd.DataFrame(bars, columns=['test', 'R2', 'P', 'n', 'colour']).drop(columns='colour').to_excel(xw, sheet_name='b_R2', index=False)
    I.to_excel(xw, sheet_name='c_indicator', index=False); dn[['sample_id', 'host_system', 'sh']].rename(columns={'sh': 'shannon'}).to_excel(xw, sheet_name='d_shannon', index=False)
    RA.to_excel(xw, sheet_name='e_relative_abundance', index=False); K.reset_index().to_excel(xw, sheet_name='f_kruskal_dunn', index=False)
print('Figure 2 written')
