"""Figure 1. Study design, sample locations and the curated gene atlas.
usage: python3 fig1_atlas.py <FINAL dir>
a workflow (SVG with host drawings; numbers from results); b world map of the 1,787 metagenomes (aggregated per host system and
0.5-degree cell; filled = site coordinates, open squares = named place or country centroid); c curated genes per family and
host system; d amino-acid identity of curated genes to the closest characterised reference.
Writes figures/Figure_1.pdf (vector), .jpg (1200 dpi), _preview.png and Figure_1_source_data.xlsx.
"""
import sys, os, json, warnings, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import *
warnings.filterwarnings('ignore')
B = sys.argv[1]; T = f'{B}/figures/tmp_fig1'; os.makedirs(T, exist_ok=True); out = f'{B}/figures/Figure_1'
S = pd.read_csv(f'{B}/results/s03_qc/samples.tsv', sep='\t'); S = S[S.set_main]
hl = json.load(open(f'{B}/results/s14_hitlevel/hitlevel_summary.json'))
mt = json.load(open(f'{B}/results/s07_metatx/metatx_summary.json'))
amb = pd.read_csv(f'{B}/results/s11_ambient/ambient_runs_rpkg.csv'); amb = amb[amb.group != 'host']
n_s, n_st = len(S), S.study.nunique()

# ---- a: workflow
icons = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'assets/host_icons.svgfrag')).read()
def box(x, y, w, h, lines, dy=34):
    t = f'<rect class="box" x="{x}" y="{y}" width="{w}" height="{h}"/>'
    y0 = y + h / 2 - (len(lines) - 1) * dy / 2 + 10
    return t + ''.join(f'<text class="bx" x="{x + w / 2}" y="{y0 + i * dy}">{l}</text>' for i, l in enumerate(lines))
svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="174mm" height="46.2mm" viewBox="0 0 1960 520" font-family="Liberation Sans, Arial, sans-serif">
<style>.hdr{{font-size:36px;font-weight:700;text-anchor:middle}} .lbl{{font-size:28px;text-anchor:middle}} .bx{{font-size:25px;fill:#111;text-anchor:middle}}
.box{{fill:#fff;stroke:#333;stroke-width:2.6;rx:7}} .panel{{fill:#fafafa;stroke:#dcdcdc;stroke-width:2;rx:14}} .arw{{stroke:#333;stroke-width:3.6;fill:none;marker-end:url(#a)}}</style>
<defs><marker id="a" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="#333"/></marker></defs>
<rect class="panel" x="20" y="92" width="750" height="350"/><rect class="panel" x="830" y="92" width="330" height="350"/>
<rect class="panel" x="1220" y="92" width="320" height="350"/><rect class="panel" x="1600" y="92" width="340" height="350"/>
<text class="hdr" x="395" y="66">Host systems</text><text class="hdr" x="995" y="66">Data</text><text class="hdr" x="1380" y="66">Gene detection</text><text class="hdr" x="1770" y="66">Analyses</text>
{icons}
<text class="lbl" x="130" y="398" fill="#E69F00">Chicken</text><text class="lbl" x="305" y="398" fill="#56B4E9">Coral</text>
<text class="lbl" x="487" y="398" fill="#009E73">Rhizosphere</text><text class="lbl" x="670" y="398" fill="#CC79A7">Sponge</text>
<path class="arw" d="M776 267 L824 267"/><path class="arw" d="M1166 267 L1214 267"/><path class="arw" d="M1546 267 L1594 267"/>
{box(848, 108, 294, 74, [f'{n_s:,} metagenomes', f'{n_st} studies'])}
{box(848, 192, 294, 74, [f"{mt['n_runs']} metatranscriptomes", f"{sum(mt['studies_by_host'].values())} studies"])}
{box(848, 276, 294, 74, [f'{amb.run.nunique()} metagenomes of', 'surroundings'])}
{box(848, 360, 294, 74, [f"{hl['genomes_searched']:,} genomes", 'from host samples'])}
{box(1236, 110, 288, 74, ['10 profile HMMs,', '217 references'])}
{box(1236, 228, 288, 74, ['phylogenetic', 'curation'])}
{box(1236, 346, 288, 80, [f"{hl['curated_genes']:,} curated genes", f"in {hl['carrier_genomes']:,} genomes"])}
<path class="arw" d="M1380 186 L1380 224"/><path class="arw" d="M1380 304 L1380 342"/>
{box(1616, 108, 308, 74, ['abundance (DNA),', 'expression (RNA)'])}
{box(1616, 192, 308, 74, ['host system and', 'study structure'])}
{box(1616, 276, 308, 74, ['subfamilies,', 'host vs surroundings'])}
{box(1616, 360, 308, 74, ['domains, structures,', 'genomic context'])}
</svg>'''
open(f'{T}/a.svg', 'w').write(svg); svg_to_pdf(svg, f'{T}/a.pdf', 174, 46.2)

# ---- b: map
import geopandas as gp
P = pd.read_csv(f'{B}/results/s08_map/sample_coordinates.tsv', sep='\t')
assert len(P) == n_s and set(P.sample_id) == set(S.sample_id), 'map table must hold exactly the main-set samples'
world = gp.read_file(gp.datasets.get_path('naturalearth_lowres')); world = world[world.name != 'Antarctica']
P['exact'] = P.coord_source.isin(['ENA run record', 'BioSample attribute'])
x = P[P.lat.notna()].copy(); x['la'] = (x.lat * 2).round() / 2; x['lo'] = (x.lon * 2).round() / 2
A = x.groupby(['host_system', 'exact', 'la', 'lo']).agg(n=('sample_id', 'size'), lat=('lat', 'mean'), lon=('lon', 'mean')).reset_index()
fig = plt.figure(figsize=(174 * MM, 80 * MM)); ax = fig.add_axes([0.0, 0.11, 1.0, 0.89])
OCEAN, LAND = '#DCE8F0', '#F3EFE4'
fig.patch.set_facecolor(OCEAN); ax.set_facecolor(OCEAN)
world.plot(ax=ax, color=LAND, edgecolor='#B9B3A3', linewidth=0.25)
size = lambda n: 3 + 6 * np.sqrt(n)
A2 = x.groupby(['host_system', 'la', 'lo']).agg(n=('sample_id', 'size'), lat=('lat', 'mean'), lon=('lon', 'mean')).reset_index()
for h in HOSTS:
    b = A2[A2.host_system == h]
    ax.scatter(b.lon, b.lat, s=size(b.n), c=COL[h], alpha=0.9, edgecolors='#333333', linewidths=0.25, zorder=3)
ax.set_xlim(-170, 180); ax.set_ylim(-66, 84); ax.set_xticks([]); ax.set_yticks([]); ax.spines[:].set_visible(False)
hs = [plt.scatter([], [], s=18, c=COL[h], label=f'{h} ({(P.host_system == h).sum()} samples, {S[S.host_system == h].study.nunique()} studies)') for h in HOSTS]
for n in [1, 10, 50]: hs.append(plt.scatter([], [], s=size(n), c='#BBBBBB', edgecolors='#333333', linewidths=0.25, label=f'{n} sample' + ('s' if n > 1 else '')))
fig.legend(handles=hs, loc='lower center', bbox_to_anchor=(0.5, 0.0), ncol=4, frameon=False, handletextpad=0.3, labelspacing=0.4, columnspacing=1.5, fontsize=6)
fig.savefig(f'{T}/b.pdf', facecolor=OCEAN); plt.close(fig)

# ---- c: curated genes per family and host
G = pd.read_csv(f'{B}/results/s14_hitlevel/genes_per_family_host.csv', index_col=0).reindex(index=FAM, columns=HOSTS)
fig, ax = plt.subplots(figsize=(104 * MM, 62 * MM)); fig.subplots_adjust(left=0.11, right=0.99, bottom=0.2, top=0.86)
w = 0.2
for j, h in enumerate(HOSTS):
    v = G[h].values; xs = np.arange(len(FAM)) + (j - 1.5) * w
    ax.bar(xs, np.where(v > 0, v, np.nan), width=w, color=COL[h], label=h)
    for xi, vi in zip(xs, v): ax.text(xi, max(vi, 1) * 1.25, f'{vi:,}' if vi > 0 else '0', ha='center', va='bottom', fontsize=4.6, rotation=90)
ax.set_yscale('log'); ax.set_ylim(0.8, 3e4); ax.set_xticks(range(len(FAM))); ax.set_xticklabels(FAM, rotation=40, ha='right', rotation_mode='anchor')
ax.set_ylabel('Curated genes'); ax.spines[['top', 'right']].set_visible(False)
ax.legend(ncol=4, frameon=False, loc='lower center', bbox_to_anchor=(0.5, 1.0), handlelength=1, columnspacing=1)
gh = pd.read_csv(f'{B}/results/s14_hitlevel/prevalence_copy_number.csv')
fig.savefig(f'{T}/c.pdf'); plt.close(fig)

# ---- d: identity to characterised references
I = pd.read_csv(f'{B}/inputs/hitlevel/identity_to_reference_per_hit.csv')
main = [f for f in FAM if f not in ('GSTDase', 'MetDase')]; hom = ['MetDase', 'GSTDase']
order = sorted(main, key=lambda f: -I[I.family == f].identity.median()) + hom
fig, ax = plt.subplots(figsize=(66 * MM, 62 * MM)); fig.subplots_adjust(left=0.36, right=0.97, bottom=0.17, top=0.86)
pos = list(range(len(main))) + [len(main) + 0.8, len(main) + 1.8]
rng = np.random.default_rng(3)
cm = plt.get_cmap('viridis'); med = I.groupby('family').identity.median()
DCOL = {f: cm(0.1 + 0.75 * (med[f] - 15) / 85) for f in main}; DCOL['MetDase'] = '#D55E00'; DCOL['GSTDase'] = '#B2182B'
for p_, f in zip(pos, order):
    v = I[I.family == f].identity.values
    cc = DCOL[f]
    bp = ax.boxplot(v, positions=[p_], vert=False, widths=0.55, showfliers=False, patch_artist=True, medianprops=dict(color='black', lw=0.9),
                    boxprops=dict(facecolor='white', edgecolor=cc, lw=0.8), whiskerprops=dict(lw=0.6, color=cc), capprops=dict(lw=0.6, color=cc))
    ax.scatter(v, p_ + rng.uniform(-0.28, 0.28, len(v)), s=1.2, c=[cc], alpha=0.35, linewidths=0, rasterized=True, zorder=1)
    for b_ in bp['boxes'] + bp['medians'] + bp['whiskers'] + bp['caps']: b_.set_zorder(4)
ax.axvline(50, color='#888888', lw=0.6, ls='--')
n_tot = I.groupby('family').size(); hits = pd.read_csv(f'{B}/results/s14_hitlevel/identity_summary.csv', index_col=0)
ax.set_yticks(pos); ax.set_yticklabels([f'{f} (n = {n_tot[f]:,})' for f in order], fontsize=5.6)
ax.invert_yaxis(); ax.set_xlim(15, 100); ax.set_xlabel('Identity to closest characterised\nreference (%)'); ax.spines[['top', 'right']].set_visible(False)
pr = hl['identity_pooled_excl_GSTDase_MetDase']
fig.savefig(f'{T}/d.pdf'); plt.close(fig)

H = assemble([('a', f'{T}/a.pdf', 0, 0, 174, 46.2), ('b', f'{T}/b.pdf', 0, 49, 174, 80),
              ('c', f'{T}/c.pdf', 0, 132, 104, 62), ('d', f'{T}/d.pdf', 106, 132, 66, 62)], out)
print('Figure 1 height mm', round(H, 1))
with pd.ExcelWriter(out + '_source_data.xlsx') as xw:
    pd.DataFrame({'item': ['metagenome samples', 'studies', 'metatranscriptome runs', 'metatranscriptome studies', 'surrounding metagenomes',
                           'genomes searched', 'curated genes', 'carrier genomes'],
                  'value': [n_s, n_st, mt['n_runs'], sum(mt['studies_by_host'].values()), amb.run.nunique(), hl['genomes_searched'], hl['curated_genes'], hl['carrier_genomes']]}
                 ).to_excel(xw, sheet_name='a_workflow', index=False)
    P.drop(columns=[c for c in ['exact', 'coord_source'] if c in P.columns]).to_excel(xw, sheet_name='b_sample_locations', index=False); A2.to_excel(xw, sheet_name='b_aggregated_points', index=False)
    G.to_excel(xw, sheet_name='c_curated_genes'); gh.to_excel(xw, sheet_name='c_prevalence_copies', index=False)
    I.to_excel(xw, sheet_name='d_identity_per_gene', index=False); hits.to_excel(xw, sheet_name='d_identity_summary')
