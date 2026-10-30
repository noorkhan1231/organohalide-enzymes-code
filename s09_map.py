"""Step 13b. World map of all samples (Figure 1 panel) and samples per country.
usage: python3 s09_map.py <sample_coordinates.tsv> <out_prefix>
Points are aggregated per host system and location (rounded to 0.5 degree); size = number of samples.
Filled circles = site coordinates (ENA run record or BioSample attribute); open squares = approximate position
(named place or country centroid). Output: <out_prefix>.pdf (vector) and <out_prefix>.jpg (1200 dpi).
"""
import sys, warnings, numpy as np, pandas as pd, geopandas as gp
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')
f, out = sys.argv[1:3]
plt.rcParams.update({'font.family': 'Liberation Sans', 'font.size': 7, 'pdf.fonttype': 42, 'svg.fonttype': 'none'})
P = pd.read_csv(f, sep='\t')
COL = {'Chicken': '#E69F00', 'Coral': '#56B4E9', 'Rhizosphere': '#009E73', 'Sponge': '#CC79A7'}   # same host colours as all manuscript figures
H = list(COL)
world = gp.read_file(gp.datasets.get_path('naturalearth_lowres')); world = world[world.name != 'Antarctica']
P['exact'] = P.coord_source.isin(['ENA run record', 'BioSample attribute'])
x = P[P.lat.notna()].copy(); x['la'] = (x.lat * 2).round() / 2; x['lo'] = (x.lon * 2).round() / 2
A = x.groupby(['host_system', 'exact', 'la', 'lo']).agg(n=('sample_id', 'size'), lat=('lat', 'mean'), lon=('lon', 'mean')).reset_index()
fig = plt.figure(figsize=(7.2, 5.6))
ax = fig.add_axes([0.0, 0.43, 1.0, 0.56])
world.plot(ax=ax, color='#E6E6E6', edgecolor='white', linewidth=0.3)
size = lambda n: 6 + 9 * np.sqrt(n)
for h in H:
    for ex in [True, False]:
        b = A[(A.host_system == h) & (A.exact == ex)]
        if not len(b): continue
        if ex: ax.scatter(b.lon, b.lat, s=size(b.n), c=COL[h], alpha=0.75, edgecolors='white', linewidths=0.3, zorder=3)
        else: ax.scatter(b.lon, b.lat, s=size(b.n), marker='s', facecolors='none', edgecolors=COL[h], linewidths=0.8, zorder=4)
ax.set_xlim(-170, 180); ax.set_ylim(-70, 84); ax.axis('off')
ax.text(0.0, 1.0, 'a', transform=ax.transAxes, fontsize=9, fontweight='bold', va='top')
nloc = len(P) - int(P.lat.notna().sum())
hs = [plt.scatter([], [], s=30, c=COL[h], label=f'{h} ({(P.host_system == h).sum()} samples, {P[P.host_system == h].study.nunique()} studies)') for h in H]
hs += [plt.scatter([], [], s=30, c='#555555', label=f'site coordinates ({int(P.exact.sum())} samples)'),
       plt.scatter([], [], s=30, marker='s', facecolors='none', edgecolors='#555555', label=f'approximate position: named place or country ({int((~P.exact & P.lat.notna()).sum())})')]
for n in [1, 10, 50]: hs.append(plt.scatter([], [], s=size(n), c='#BBBBBB', label=f'{n} sample' + ('s' if n > 1 else '')))
leg = fig.legend(handles=hs, loc='upper center', bbox_to_anchor=(0.5, 0.435), fontsize=6, frameon=False, ncol=3, handletextpad=0.4, labelspacing=0.5, columnspacing=1.2)
# panel b: samples per country / sea region (top 25)
c = P.country.fillna('').replace('', 'not recorded')
c = c.replace({'United States of America': 'USA', 'Viet Nam': 'Vietnam', 'United Arab Emirates': 'UAE', 'United Kingdom': 'UK'})
T = pd.crosstab(c, P.host_system).reindex(columns=H, fill_value=0)
T = T.loc[T.sum(axis=1).sort_values(ascending=False).index]
T = T[T.index != 'not recorded'].head(25)
bx = fig.add_axes([0.07, 0.11, 0.92, 0.19])
bottom = np.zeros(len(T))
for h in H:
    bx.bar(range(len(T)), T[h].values, bottom=bottom, color=COL[h], width=0.75, linewidth=0); bottom += T[h].values
bx.set_xticks(range(len(T))); bx.set_xticklabels(T.index, rotation=45, ha='right', rotation_mode='anchor', fontsize=6)
bx.set_ylabel('Samples'); bx.spines[['top', 'right']].set_visible(False); bx.set_xlim(-0.6, len(T) - 0.4)
bx.text(-0.06, 1.05, 'b', transform=bx.transAxes, fontsize=9, fontweight='bold', va='bottom')
fig.savefig(out + '.pdf'); fig.savefig(out + '.jpg', dpi=1200, pil_kwargs={'quality': 95}); fig.savefig(out + '_preview.png', dpi=200)
T.to_csv(out + '_samples_per_country.csv')
print('countries/sea regions with samples:', int((c != 'not recorded').groupby(c).size().shape[0] - (1 if 'not recorded' in set(c) else 0)), '| location not recorded:', nloc)
