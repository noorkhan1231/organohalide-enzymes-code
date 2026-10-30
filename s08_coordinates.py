"""Step 13a. Coordinates for every sample in the main set.
usage: python3 s08_coordinates.py <inputs_dir> <s03 dir> <out_dir>
Order of sources: (1) ENA run record lat/lon; (2) BioSample attributes lat_lon or geographic location (latitude/longitude)
from NCBI or ENA sample records; (3) country (or sea) centroid from the ENA country field, flagged 'approximate'.
Values outside [-90, 90] / [-180, 180] or written as missing/not applicable are ignored.
"""
import sys, os, re, glob, pandas as pd, numpy as np
I, Q, O = sys.argv[1:4]; os.makedirs(O, exist_ok=True)
R = pd.read_csv(f'{Q}/run_qc.tsv', sep='\t', low_memory=False); S = pd.read_csv(f'{Q}/samples.tsv', sep='\t'); S = S[S.set_main]
inc = R[R.sample_id.isin(S.sample_id)]
def num(x):
    try: v = float(x); return v if np.isfinite(v) else np.nan
    except Exception: return np.nan
def parse_latlon(v):
    v = str(v).strip()
    m = re.match(r'^\s*(-?\d+(?:\.\d+)?)\s*([NS])?[\s,;]+(-?\d+(?:\.\d+)?)\s*([EW])?\s*$', v, re.I)
    if not m: return np.nan, np.nan
    la, lo = float(m.group(1)), float(m.group(3))
    if m.group(2) and m.group(2).upper() == 'S': la = -abs(la)
    if m.group(4) and m.group(4).upper() == 'W': lo = -abs(lo)
    return la, lo
# BioSample attribute files
bs = []
for f in [f'{I}/metadata/biosample_geo_original.tsv', f'{I}/metadata/biosample_attributes_new.tsv', f'{I}/metadata/bs_extra_results.tsv']:
    if os.path.exists(f): bs.append(pd.read_csv(f, sep='\t', dtype=str))
B = pd.concat(bs, ignore_index=True).dropna(subset=['attribute', 'value'])
bll = {}
for s, g in B.groupby('sample_accession'):
    a = {k.lower(): v for k, v in zip(g.attribute, g.value)}
    la = lo = np.nan
    for k in ['lat_lon', 'lat lon', 'geographic location (latitude and longitude)']:
        if k in a: la, lo = parse_latlon(a[k]); break
    if np.isnan(la):
        la = num(re.sub(r'[^\d\.\-]', '', a.get('geographic location (latitude)', '').split(' ')[0])) if a.get('geographic location (latitude)') else np.nan
        lo = num(re.sub(r'[^\d\.\-]', '', a.get('geographic location (longitude)', '').split(' ')[0])) if a.get('geographic location (longitude)') else np.nan
        if 'S' in a.get('geographic location (latitude)', '') and la > 0: la = -la
        if 'W' in a.get('geographic location (longitude)', '') and lo > 0: lo = -lo
    if -90 <= la <= 90 and -180 <= lo <= 180 and not (la == 0 and lo == 0): bll[s] = (la, lo)
rows = []
for _, r in inc.drop_duplicates('sample_id').iterrows():
    la, lo, src = num(r.lat), num(r.lon), 'ENA run record'
    if not (-90 <= la <= 90 and -180 <= lo <= 180) or (la == 0 and lo == 0):
        la = lo = np.nan
        if r.sample_accession in bll: (la, lo), src = bll[r.sample_accession], 'BioSample attribute'
    rows.append(dict(sample_id=r.sample_id, host_system=r.host_system, BioProject=r.BioProject, study=r.study, country_field=r.country, lat=la, lon=lo, coord_source=src if pd.notna(la) else ''))
P = pd.DataFrame(rows)
P['country'] = P.country_field.fillna('').astype(str).str.split(':').str[0].str.strip()
P.loc[P.country.isin(['', 'missing', 'not applicable', 'not collected', 'nan']), 'country'] = ''
P.to_csv(f'{O}/sample_coordinates_raw.tsv', sep='\t', index=False)
print('samples', len(P), '| with site coordinates', int(P.lat.notna().sum()), P[P.lat.notna()].coord_source.value_counts().to_dict())
print('without coordinates by host:', P[P.lat.isna()].host_system.value_counts().to_dict())
print('without coordinates: country field values:'); print(P[P.lat.isna()].groupby(['host_system', 'country']).size().to_string())

# (4) named places (manual, approximate) and (5) country centroids (approximate)
import geopandas as gp
M = pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'place_coordinates_manual.tsv'), sep='\t')
mp = M.set_index('match_country_field')
k = P.lat.isna() & P.country_field.isin(mp.index)
P.loc[k, 'lat'] = P.loc[k, 'country_field'].map(mp.lat); P.loc[k, 'lon'] = P.loc[k, 'country_field'].map(mp.lon); P.loc[k, 'coord_source'] = 'approximate: named place'
world = gp.read_file(gp.datasets.get_path('naturalearth_lowres'))
NE = {'USA': 'United States of America', 'United States': 'United States of America', 'UK': 'United Kingdom', 'Czech Republic': 'Czechia', 'Viet Nam': 'Vietnam', 'Russia': 'Russia', 'South Korea': 'South Korea'}
def _c(g):
    c = g.centroid
    if g.contains(c): return c
    big = max(g.geoms, key=lambda p: p.area) if g.geom_type == 'MultiPolygon' else g     # centroid of the main landmass
    c = big.centroid
    return c if big.contains(c) else big.representative_point()
cent = {r.name: _c(r.geometry) for r in world.itertuples()}
FIX = {'PRJEB38174': 'France'}       # MetaChick: ENA country empty; BioSample geo_loc_name France (Evry = sequencing centre, not farm) -> country centroid
P.loc[P.lat.isna() & P.BioProject.isin(FIX), 'country'] = P.BioProject.map(FIX)
k = P.lat.isna() & P.country.map(lambda c: NE.get(c, c) in cent)
P.loc[k, 'lat'] = P.loc[k, 'country'].map(lambda c: cent[NE.get(c, c)].y); P.loc[k, 'lon'] = P.loc[k, 'country'].map(lambda c: cent[NE.get(c, c)].x)
P.loc[k, 'coord_source'] = 'approximate: country centroid'
P.loc[P.lat.isna(), 'coord_source'] = 'location not recorded'
P.to_csv(f'{O}/sample_coordinates.tsv', sep='\t', index=False)
print('\nfinal coordinate sources:'); print(P.groupby(['host_system', 'coord_source']).size().unstack(fill_value=0).to_string())
