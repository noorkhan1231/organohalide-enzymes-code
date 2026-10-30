"""Step 6. Sample and study tables for the manuscript (Supplementary Table: samples and BioProjects per environment).
usage: python3 s03b_tables.py <s03 dir> <out.xlsx>"""
import sys, pandas as pd, numpy as np
Q, out = sys.argv[1:3]
R = pd.read_csv(f'{Q}/run_qc.tsv', sep='\t', low_memory=False); S = pd.read_csv(f'{Q}/samples.tsv', sep='\t'); S = S[S.set_main]
H = ['Chicken', 'Coral', 'Rhizosphere', 'Sponge']
inc = R[R.sample_id.isin(S.sample_id)].copy()
cty = lambda s: '; '.join(sorted({str(x).split(':')[0].strip() for x in s.dropna() if str(x).strip() and str(x) not in ('missing', 'not applicable', 'nan')})) or 'not recorded'
top = lambda s: '; '.join(f'{k} ({v})' for k, v in s.dropna().astype(str).str.strip().replace('', np.nan).dropna().value_counts().head(3).items())
host = inc.host.where(inc.host.notna() & ~inc.host.astype(str).str.lower().isin(['missing', 'not applicable', 'na', 'not provided']), inc.scientific_name)
HOST_FIX = {'PRJNA1231735': 'Glycine max (rhizoplane; BioSample env_local_scale soybean)', 'PRJNA688820': 'Triticum aestivum (6 cultivars; BioSample cultivar field)',
            'PRJEB52368': 'Tara Pacific coral holobiont (BioSample environment C-CORAL)'}
host = host.where(~inc.BioProject.isin(HOST_FIX), inc.BioProject.map(HOST_FIX))
inc['host_or_organism'] = host
B = inc.groupby(['host_system', 'study', 'BioProject']).agg(study_title=('study_title', 'first'), runs=('run', 'size'), samples=('sample_id', 'nunique'),
        country=('country', cty), host_species_or_organism=('host_or_organism', top), sample_type=('isolation_source', top),
        median_Gb=('Gb', 'median')).reset_index()
B = B.merge(S.groupby('BioProject').GE.median().rename('median_genome_equivalents').reset_index(), on='BioProject', how='left')
B = B.sort_values(['host_system', 'samples'], ascending=[True, False])
Sm = pd.DataFrame([dict(environment=h, samples=int((S.host_system == h).sum()), runs=int((inc.host_system == h).sum()),
        BioProjects=int(S[S.host_system == h].BioProject.nunique()), independent_studies=int(S[S.host_system == h].study.nunique()),
        countries=int(len({c for c in '; '.join(B[B.host_system == h].country).split('; ') if c != 'not recorded'})),
        median_samples_per_study=float(S[S.host_system == h].groupby('study').size().median()),
        max_samples_one_study=int(S[S.host_system == h].groupby('study').size().max()),
        median_genome_equivalents=float(S[S.host_system == h].GE.median())) for h in H])
Sm.loc[len(Sm)] = dict(environment='Total', samples=len(S), runs=len(inc), BioProjects=S.BioProject.nunique(), independent_studies=S.study.nunique(),
        countries=len({c for c in '; '.join(B.country).split('; ') if c != 'not recorded'}), median_samples_per_study=float(S.groupby('study').size().median()),
        max_samples_one_study=int(S.groupby('study').size().max()), median_genome_equivalents=float(S.GE.median()))
Ss = S[['sample_id', 'host_system', 'BioProject', 'study', 'runs', 'n_runs', 'country', 'lat', 'lon', 'reads_processed', 'Gb', 'GE', 'GE_per_M', 'platform']].sort_values(['host_system', 'BioProject', 'sample_id'])
X = R[R.exclude_code.fillna('') != ''][['run', 'host_system', 'BioProject', 'sample_accession', 'exclude_code', 'exclude_reason']].sort_values(['host_system', 'exclude_code', 'BioProject'])
XS = X.groupby(['host_system', 'exclude_code']).size().unstack(fill_value=0)
with pd.ExcelWriter(out) as w:
    Sm.to_excel(w, sheet_name='Summary', index=False); B.to_excel(w, sheet_name='BioProjects', index=False); Ss.to_excel(w, sheet_name='Samples', index=False)
    XS.to_excel(w, sheet_name='Excluded_summary'); X.to_excel(w, sheet_name='Excluded_runs', index=False)
print(Sm.to_string(index=False)); print(XS.to_string())
B.to_csv(out.replace('.xlsx', '_bioprojects.tsv'), sep='\t', index=False)
