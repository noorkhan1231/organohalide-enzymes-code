"""Figure 5. Domain architecture, predicted structures and genomic context of the curated genes.
usage: python3 fig5_architecture.py <FINAL dir>
a Pfam-A architectures (schematic, not to scale) of the RDase CxxCxxCxxxCP motif group and of VanHase, short and extended forms;
b share of genes carrying the extra module per host system; c domains that differ among host systems beyond GTDB order
(share of carrier genomes with the domain); d genes on geNomad plasmid-classified contigs, observed vs expected from contig-length-matched background; e pLDDT and f best Foldseek TM-score of the 529 predicted structures.
"""
import sys, os, json, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import *
from matplotlib.patches import Rectangle, Patch
B = sys.argv[1]; T = f'{B}/figures/tmp_fig5'; os.makedirs(T, exist_ok=True); out = f'{B}/figures/Figure_5'
A = pd.read_excel(f'{B}/inputs/hitlevel/Domain_architecture_per_hit.xlsx'); A['doms'] = A.arch.fillna('').str.split('|')
DV = pd.read_excel(f'{B}/inputs/hitlevel/Domain_vs_host_all_hits.xlsx', 'domain_vs_host')
RM = pd.read_csv(f'{B}/results/s14_hitlevel/rdase_module.csv', index_col=0); VM = pd.read_csv(f'{B}/results/s14_hitlevel/vanhase_module.csv', index_col=0)
ST = pd.read_csv(f'{B}/results/s13_structures/structures_529_per_structure.csv')
MB = pd.read_csv(f'{B}/results/s14_hitlevel/mobility_plasmid_contig.csv')
HL = json.load(open(f'{B}/results/s14_hitlevel/hitlevel_summary.json'))
DT = pd.read_csv(f'{B}/inputs/hitlevel/deeptmhmm_per_hit_annotated.tsv', sep='\t').set_index('hit_id')

# ---- a: architecture schematics
DC = {'Fer4_7': '#4C72B0', 'Fer4_16': '#6A8FC7', 'Fer2': '#C44E52', 'DmmA-like_N': '#E39A9C', 'TAT': '#999999', 'SP': '#999999',
      'VCPO_2nd': '#55A868', 'PAP2': '#8BC59B', 'DUF6851': '#B5D9BC', 'FG-GAP_3': '#8172B2', 'UnbV_ASPIC': '#A99BD0', 'Flg_new_2': '#CFC6E8'}
r = A[(A.gene == 'RDase') & (A.subfamily_call == 'motifgroup_CxxCxxCxxxCP')].copy(); r['mod'] = r.doms.apply(lambda d: 'Fer2' in d or 'DmmA-like_N' in d)
v = A[A.gene == 'VanHase'].copy(); v['mod'] = v.doms.apply(lambda d: 'FG-GAP_3' in d or 'UnbV_ASPIC' in d)
ln = lambda ids: int(DT.loc[ids, 'len'].median())
SUB = [r[~r['mod']], r[r['mod']], v[~v['mod']], v[v['mod']]]
rows = [('RDase motif group', ['Fer4_7', 'Fer4_16'], ln(r[~r['mod']].hit_id), (~r['mod']).sum()),
        ('RDase motif group, extended', ['Fer4_7', 'Fer4_16', 'DmmA-like_N', 'Fer2'], ln(r[r['mod']].hit_id), r['mod'].sum()),
        ('VanHase', ['VCPO_2nd'], ln(v[~v['mod']].hit_id), (~v['mod']).sum()),
        ('VanHase, extended', ['SP', 'FG-GAP_3', 'UnbV_ASPIC', 'Flg_new_2', 'VCPO_2nd'], ln(v[v['mod']].hit_id), v['mod'].sum())]
sp_ext = int(DT.loc[v[v['mod']].hit_id, 'SP'].sum())
fig, ax = plt.subplots(figsize=(104 * MM, 50 * MM)); fig.subplots_adjust(left=0.01, right=0.99, bottom=0.08, top=0.98)
ax.set_xlim(0, 100); ax.set_ylim(4.5, -1.0); ax.axis('off')
ax.text(90, -0.75, 'Host system share', ha='center', va='center', fontsize=5.4, fontweight='bold')
for i, (name, doms, L, n) in enumerate(rows):
    ev = SUB[i].environment.value_counts().reindex(HOSTS).fillna(0); xs = 82.0
    for h in HOSTS:
        wd = 16.0 * ev[h] / ev.sum()
        if wd > 0: ax.add_patch(Rectangle((xs, i - 0.13), wd, 0.26, facecolor=COL[h], edgecolor='white', lw=0.3, zorder=2)); xs += wd
    ax.text(3, i, f'{name}\n(n = {n}; median {L} aa)', va='center', fontsize=5.6)
    x0, x1 = 34, 34 + 44 * min(L, 1300) / 1300
    ax.plot([x0, x1], [i, i], color='#555555', lw=1.0, solid_capstyle='butt', zorder=1)
    w = (x1 - x0) / len(doms)
    for k, d in enumerate(doms):
        bw = w * (0.35 if d in ('SP', 'TAT') else 0.82)
        ax.add_patch(Rectangle((x0 + k * w + (w - bw) / 2, i - 0.26), bw, 0.52, facecolor=DC[d], edgecolor='none', zorder=2))
        ax.text(x0 + k * w + w / 2, i + (0.42 if (k % 2 == 0 or len(doms) < 5) else 0.62), d, ha='center', va='top', fontsize=4.4)
ax.legend(handles=[Patch(facecolor=COL[h], label=h) for h in HOSTS], ncol=4, frameon=False, loc='lower left', bbox_to_anchor=(0.0, -0.04), handlelength=0.9, columnspacing=1.0, fontsize=5.2)
fig.savefig(f'{T}/a.pdf'); plt.close(fig)

# ---- b: module carriage per host
fig, axs = plt.subplots(1, 2, figsize=(68 * MM, 50 * MM)); fig.subplots_adjust(left=0.16, right=0.98, bottom=0.27, top=0.86, wspace=0.35)
for ax, (tab, title) in zip(axs, [(RM, 'RDase motif group\nwith Fer2/DmmA-like N'), (VM, 'VanHase with\nFG-GAP/UnbV')]):
    pct = 100 * tab.with_module / tab.genes.replace(0, np.nan)
    ax.bar(range(4), pct.fillna(0).values, color=[COL[h] for h in HOSTS], width=0.7)
    for j, h in enumerate(HOSTS): ax.text(j, (pct.fillna(0).iloc[j]) + 2, f'{tab.with_module[h]}/{tab.genes[h]}', ha='center', fontsize=4.8, rotation=90, va='bottom')
    ax.set_xticks(range(4)); ax.set_xticklabels(HOSTS, rotation=90, fontsize=5.4); ax.set_ylim(0, 100); ax.set_title(title, fontsize=5.6)
    ax.spines[['top', 'right']].set_visible(False)
axs[0].set_ylabel('Genes (%)')
fig.savefig(f'{T}/b.pdf'); plt.close(fig)

# ---- c: domains beyond lineage (genome-level presence within family)
D = DV[DV.call == 'host-associated beyond lineage'].copy()
M = []
for d in D.itertuples():
    x = A[A.gene == d.gene].copy(); x['has'] = x.doms.apply(lambda L: d.pfam in L)
    g = x.groupby(['environment', 'genome_id']).has.any().groupby('environment').agg(['mean', 'size'])
    M.append(dict(label=f'{d.gene} {d.pfam}', **{h: 100 * g['mean'].get(h, np.nan) if g['size'].get(h, 0) >= 5 else np.nan for h in HOSTS}, FDR_within_order=d.FDR_within_order))
M = pd.DataFrame(M).set_index('label')
fig, ax = plt.subplots(figsize=(88 * MM, 66 * MM)); fig.subplots_adjust(left=0.36, right=0.88, bottom=0.2, top=0.98)
im = ax.imshow(M[HOSTS].values, cmap='YlGnBu', vmin=0, vmax=100, aspect='auto')
for i in range(len(M)):
    for j in range(4):
        val = M[HOSTS].values[i, j]
        ax.text(j, i, '–' if np.isnan(val) else f'{val:.0f}', ha='center', va='center', fontsize=5.0, color='white' if (not np.isnan(val) and val > 60) else 'black')
ax.set_xticks(range(4)); ax.set_xticklabels(HOSTS, rotation=35, ha='right', rotation_mode='anchor'); ax.set_yticks(range(len(M))); ax.set_yticklabels(M.index, fontsize=5.4)
cb = fig.colorbar(im, ax=ax, fraction=0.06, pad=0.03); cb.set_label('Carrier genomes with domain (%)', fontsize=5.6); cb.ax.tick_params(labelsize=5.2)
fig.savefig(f'{T}/c.pdf'); plt.close(fig)

# ---- f: plasmid-classified contigs
fig, ax = plt.subplots(figsize=(80 * MM, 66 * MM)); fig.subplots_adjust(left=0.17, right=0.98, bottom=0.2, top=0.92)
lab = [f'{r.family}\n{r.host_system}' for r in MB.itertuples()]; xx = np.arange(len(MB))
ax.bar(xx - 0.18, MB.plasmid_pct, width=0.36, color='#0072B2', label='observed')
ax.bar(xx + 0.18, MB.expected_pct, width=0.36, color='#F0C419', label='expected')
for k, r in enumerate(MB.itertuples()):
    s = '***' if r.FDR < 0.001 else '**' if r.FDR < 0.01 else '*' if r.FDR < 0.05 else 'n.s.'
    ax.text(k, max(r.plasmid_pct, r.expected_pct) + 0.5, s, ha='center', fontsize=5.6)
ax.set_xticks(xx); ax.set_xticklabels(lab, fontsize=5.2); ax.set_ylabel('Genes on plasmid-classified contigs (%)')
ax.legend(handles=[Patch(facecolor='#0072B2', label='observed'), Patch(facecolor='#F0C419', label='expected')], frameon=False, loc='upper left', handlelength=1)
ax.spines[['top', 'right']].set_visible(False); ax.set_ylim(0, 18)
fig.savefig(f'{T}/f.pdf'); plt.close(fig)

# ---- d, e: structures
fams = [f for f in FAM if f in set(ST.gene)]
order = sorted(fams, key=lambda f: -ST[ST.gene == f].shape[0])
rng = np.random.default_rng(11)
for key, col, lim, ref, ylab in [('d', 'plddt', (60, 100), 70, 'Mean pLDDT'), ('e', 'best_tm_score', (0.4, 1.02), 0.5, 'Best TM-score (Foldseek)')]:
    fig, ax = plt.subplots(figsize=(86 * MM, 48 * MM)); fig.subplots_adjust(left=0.14, right=0.98, bottom=0.22, top=0.9)
    for i, f in enumerate(order):
        x = ST[ST.gene == f]
        ax.scatter(i + rng.uniform(-0.28, 0.28, len(x)), x[col], s=2.2, c=[COL[h] for h in x.environment], linewidths=0, alpha=0.8, rasterized=True)
        ax.plot([i - 0.32, i + 0.32], [x[col].median()] * 2, color='black', lw=0.9)
    ax.axhline(ref, color='#888888', lw=0.6, ls='--'); ax.set_ylim(*lim)
    ax.set_xticks(range(len(order))); ax.set_xticklabels([f'{f}\n({(ST.gene == f).sum()})' for f in order], fontsize=5.2)
    ax.set_ylabel(ylab); ax.spines[['top', 'right']].set_visible(False)
    if True:
        ax.legend(handles=[Patch(facecolor=COL[h], label=h) for h in HOSTS], ncol=4, frameon=False, loc='lower center', bbox_to_anchor=(0.5, 0.98), handlelength=0.8, columnspacing=0.8, fontsize=5.4)
    fig.savefig(f'{T}/{key}.pdf'); plt.close(fig)

H = assemble([('a', f'{T}/a.pdf', 0, 0, 104, 50), ('b', f'{T}/b.pdf', 106, 0, 68, 50),
              ('c', f'{T}/c.pdf', 0, 53, 88, 66), ('d', f'{T}/f.pdf', 92, 53, 80, 66),
              ('e', f'{T}/d.pdf', 0, 122, 86, 48), ('f', f'{T}/e.pdf', 88, 122, 86, 48)], out)
print('Figure 5 height mm', round(H, 1), '| extended VanHase with signal peptide', sp_ext, 'of', v['mod'].sum())
with pd.ExcelWriter(out + '_source_data.xlsx') as xw:
    pd.DataFrame([(n_, '|'.join(d_), L_, k_) for n_, d_, L_, k_ in rows], columns=['architecture', 'domains_N_to_C_schematic', 'median_length_aa', 'genes']).to_excel(xw, sheet_name='a_architectures', index=False)
    RM.to_excel(xw, sheet_name='b_rdase_module'); VM.to_excel(xw, sheet_name='b_vanhase_module')
    M.to_excel(xw, sheet_name='c_domain_presence_pct'); DV.to_excel(xw, sheet_name='c_all_domain_tests', index=False)
    ST.to_excel(xw, sheet_name='ef_structures_529', index=False); MB.to_excel(xw, sheet_name='d_plasmid_contigs', index=False)
