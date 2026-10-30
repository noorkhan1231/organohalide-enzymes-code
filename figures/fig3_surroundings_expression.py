"""Figure 3. Host-associated communities versus their surroundings, and expression of the ten families.
usage: python3 fig3_surroundings_expression.py <FINAL dir>
a reads recruited by random non-target (null) genes of each system's genome collection, host vs surroundings (study means);
b host vs surroundings call per family (agreement of RPKG and null-normalised abundance, both FDR < 0.05; number = null-normalised
  host/surroundings ratio of study means); c share of metatranscriptome runs with reads on each family; d expression (reads per
  kb per million reads processed) per run, with run- and study-level Kruskal-Wallis FDR.
"""
import sys, os, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import *
from matplotlib.patches import Rectangle, Patch
B = sys.argv[1]; T = f'{B}/figures/tmp_fig3'; os.makedirs(T, exist_ok=True); out = f'{B}/figures/Figure_3'
R = pd.read_csv(f'{B}/results/s11_ambient/ambient_runs_rpkg.csv')
C = pd.read_csv(f'{B}/results/s11_ambient/host_vs_ambient_calls.csv')
NL = pd.read_csv(f'{B}/results/s11_ambient/null_gene_recruitment.csv').set_index('system')
W = pd.read_csv(f'{B}/results/s07_metatx/metatx_rpkm_per_run.csv'); K = pd.read_csv(f'{B}/results/s07_metatx/metatx_family_tests.csv').set_index('family')
GREY = '#D55E00'; NOCOL = '#F6F0C4'; UNCOL = '#C9B8E3'

# ---- a: null-gene recruitment, study means
x = R[R.null_RPKG > 0].copy(); x['side'] = np.where(x.group == 'host', 'host', 'surroundings')
st = x.groupby(['system', 'side', 'unit']).null_RPKG.mean().reset_index()
fig = plt.figure(figsize=(70 * MM, 72 * MM)); ax = fig.add_axes([0.17, 16 / 72, 0.81, 42 / 72])
rng = np.random.default_rng(5)
for i, h in enumerate(HOSTS):
    for k, side in enumerate(['host', 'surroundings']):
        v = st[(st.system == h) & (st.side == side)].null_RPKG.values; p = i + (k - 0.5) * 0.38
        c = COL[h] if side == 'host' else GREY
        ax.boxplot(np.log10(v), positions=[p], widths=0.3, showfliers=False, patch_artist=True, medianprops=dict(color='black', lw=0.7),
                   boxprops=dict(facecolor=c, alpha=0.55, lw=0.5), whiskerprops=dict(lw=0.5), capprops=dict(lw=0.5))
        ax.scatter(p + rng.uniform(-0.1, 0.1, len(v)), np.log10(v), s=2, c=c, linewidths=0, zorder=3)
    r = NL.loc[h]; ax.text(i, 3.45, f"{r.ratio:.1f}×\n" + ('FDR <0.001' if r.FDR < 0.001 else f'FDR {r.FDR:.3f}'), ha='center', va='bottom', fontsize=4.8)
ax.set_xticks(range(4)); ax.set_xticklabels([f'{h}\n{int(NL.loc[h, "n_host_studies"])} vs {int(NL.loc[h, "n_ambient_studies"])}' for h in HOSTS], fontsize=5.6)
ax.set_ylim(-3.2, 4.6); ax.set_yticks(range(-3, 5)); ax.set_yticklabels([f'$10^{{{k}}}$' for k in range(-3, 5)])
ax.set_ylabel('Null-gene RPKG (study mean)'); ax.spines[['top', 'right']].set_visible(False)
fig.savefig(f'{T}/a.pdf'); plt.close(fig)

# ---- b: call matrix
fig = plt.figure(figsize=(102 * MM, 72 * MM)); ax = fig.add_axes([0.17, 16 / 72, 0.82, 42 / 72])
def ratio_txt(r):
    if np.isinf(r): return 'host only'
    if np.isnan(r): return '–'
    if r < 0.01: return '<0.01×'
    return f'{r:.2f}×' if r < 10 else f'{r:.0f}×'
for i, f in enumerate(FAM):
    for j, h in enumerate(HOSTS):
        r = C[(C.system == h) & (C.set == f)].iloc[0]; cl = r.call
        fc, tc, hatch = NOCOL, '#555555', None
        if cl == 'host-enriched': fc, tc = COL[h], 'black'
        elif cl == 'surroundings-enriched': fc, tc = GREY, 'black'
        elif cl == 'unresolved': fc, tc, hatch = UNCOL, '#222222', '////'
        matplotlib.rcParams['hatch.linewidth'] = 0.4
        ax.add_patch(Rectangle((j, i), 1, 1, facecolor=fc, edgecolor='white', lw=0.8, hatch=hatch))
        if r.host_studies_detected == 0 and r.ambient_studies_detected == 0: t = 'not detected'
        else: t = ratio_txt(r.ratio_null_normalised)
        ax.text(j + 0.5, i + 0.5, t, ha='center', va='center', fontsize=5.3, color=tc)
ax.set_xlim(0, 4); ax.set_ylim(len(FAM), 0); ax.set_yticks(np.arange(len(FAM)) + 0.5); ax.set_yticklabels(FAM)
c0 = C.groupby('system')[['n_host_studies', 'n_ambient_studies']].max()
ax.set_xticks(np.arange(4) + 0.5); ax.xaxis.tick_top()
ax.set_xticklabels([f'{h}\n{c0.loc[h, "n_host_studies"]} vs {c0.loc[h, "n_ambient_studies"]}' for h in HOSTS], fontsize=5.8)
ax.tick_params(length=0); [s_.set_visible(False) for s_ in ax.spines.values()]
lg = [Patch(facecolor='#56B4E9', edgecolor='white', label='Host-enriched (host-system colour)'), Patch(facecolor=GREY, label='Surroundings-enriched'),
      Patch(facecolor=UNCOL, edgecolor='#7B5EA7', hatch='////', label='Unresolved (RPKG and null-normalised abundance disagree)'), Patch(facecolor=NOCOL, label='No consistent difference')]
fig.legend(handles=lg, loc='lower left', bbox_to_anchor=(0.17, 0.0), frameon=False, fontsize=5.6, handlelength=1.4, labelspacing=0.5, ncol=1)
fig.savefig(f'{T}/b.pdf'); plt.close(fig)

# ---- c: detection in metatranscriptomes
D = pd.DataFrame({h: [(W.loc[W.host_system == h, f] > 0).mean() * 100 for f in FAM] for h in HOSTS}, index=FAM)
nrun = W.host_system.value_counts(); nst = W.groupby('host_system').study.nunique()
fig, ax = plt.subplots(figsize=(56 * MM, 66 * MM)); fig.subplots_adjust(left=0.25, right=0.84, bottom=0.27, top=0.98)
im = ax.imshow(D.values, cmap='Purples', vmin=0, vmax=100, aspect='auto')
for i in range(len(FAM)):
    for j in range(4):
        v = D.values[i, j]; ax.text(j, i, f'{v:.0f}', ha='center', va='center', fontsize=5.4, color='white' if v > 55 else 'black')
ax.set_xticks(range(4)); ax.set_xticklabels([f'{h} ({nrun[h]}, {nst[h]})' for h in HOSTS], fontsize=5.4, rotation=40, ha='right', rotation_mode='anchor')
ax.set_yticks(range(len(FAM))); ax.set_yticklabels(FAM)
assert (K.FDR_detection < 0.05).all()
cb = fig.colorbar(im, ax=ax, fraction=0.08, pad=0.03); cb.set_label('Runs (%)', fontsize=5.6, labelpad=2); cb.ax.tick_params(labelsize=5.2)
fig.savefig(f'{T}/c.pdf'); plt.close(fig)

# ---- d: expression level
fig, ax = plt.subplots(figsize=(116 * MM, 66 * MM)); fig.subplots_adjust(left=0.13, right=0.99, bottom=0.24, top=0.86)
w = 0.2; floor = 1e-3
for j, h in enumerate(HOSTS):
    for i, f in enumerate(FAM):
        v = W.loc[W.host_system == h, f].values; vp = v[v > 0]; p = i + (j - 1.5) * w
        if len(vp) >= 3:
            ax.boxplot(np.log10(vp), positions=[p], widths=0.17, showfliers=False, patch_artist=True, medianprops=dict(color='black', lw=0.6),
                       boxprops=dict(facecolor=COL[h], alpha=0.6, lw=0.4), whiskerprops=dict(lw=0.4), capprops=dict(lw=0.4))
        if len(vp): ax.scatter(p + rng.uniform(-0.06, 0.06, len(vp)), np.log10(vp), s=1.2, c=COL[h], linewidths=0, zorder=3, rasterized=True)
for i, f in enumerate(FAM):
    k = K.loc[f]
    s1 = '*' if k.FDR < 0.05 else 'ns'; s2 = '†' if k.FDR_study < 0.05 else ''
    ax.text(i, 3.25, s1 + s2, ha='center', fontsize=6.5)
ax.set_xticks(range(len(FAM))); ax.set_xticklabels(FAM, rotation=40, ha='right', rotation_mode='anchor')
ax.set_ylim(-2.2, 3.6); ax.set_yticks(range(-2, 4)); ax.set_yticklabels([f'$10^{{{k}}}$' for k in range(-2, 4)])
ax.set_ylabel('Reads per kb per\nmillion reads'); ax.spines[['top', 'right']].set_visible(False)
ax.legend(handles=[Patch(facecolor=COL[h], alpha=0.7, label=h) for h in HOSTS], ncol=4, frameon=False, loc='lower center', bbox_to_anchor=(0.5, 1.02), handlelength=1)
fig.savefig(f'{T}/d.pdf'); plt.close(fig)

H = assemble([('a', f'{T}/a.pdf', 0, 0, 70, 72), ('b', f'{T}/b.pdf', 72, 0, 102, 72),
              ('c', f'{T}/c.pdf', 0, 75, 56, 66), ('d', f'{T}/d.pdf', 58, 75, 116, 66)], out)
print('Figure 3 height mm', round(H, 1))
with pd.ExcelWriter(out + '_source_data.xlsx') as xw:
    st.to_excel(xw, sheet_name='a_null_study_means', index=False); NL.reset_index().to_excel(xw, sheet_name='a_null_test', index=False)
    C.to_excel(xw, sheet_name='b_calls', index=False)
    pd.read_csv(f'{B}/results/s11_ambient/host_vs_ambient_family.csv').to_excel(xw, sheet_name='b_all_tests', index=False)
    D.to_excel(xw, sheet_name='c_detection_pct'); W.to_excel(xw, sheet_name='d_rpkm_per_run', index=False); K.reset_index().to_excel(xw, sheet_name='cd_tests', index=False)
