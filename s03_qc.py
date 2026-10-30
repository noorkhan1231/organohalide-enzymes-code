"""Step 5-6. Sample QC for every metagenome run, BioSample pooling, genome-equivalent floor and study groups.
usage: python3 s03_qc.py <s02 metadata table> <qc_decisions.tsv> <study_groups.tsv> <out_dir> [GE_floor=50]
Rules, applied identically to all runs (original and added):
  LIB     ENA library_strategy is not WGS (amplicon, WGA, Hi-C, OTHER, pool-clone libraries)
  NOENA   run absent from ENA (metadata cannot be verified)
  decision table (qc_decisions.tsv): ENV environmental sample (seawater, sediment, bulk soil); FW freshwater sponge
          (decision D1); VIR virus-enriched fraction; HOSTGEN host-genome / genome-skimming study; OTHERHOST wrong host
          organism; TISSUE compartment outside scope (root endosphere, seed-borne, FFPE micro-dissected tissue);
          CAPTIVE home aquarium without provenance; DUP possible re-sequencing of specimens already included
  pooling runs of one BioSample are pooled (reads, genome equivalents and target counts summed) -> one sample
  LOWGE   pooled genome equivalents < floor (main set); sensitivity sets keep GE >= 20 and GE >= 200
Outputs: run_qc.tsv (every run + decision), samples.tsv (one row per pooled sample in the analysis sets), summary.
"""
import sys, os, pandas as pd, numpy as np
mf, df, gf, O = sys.argv[1:5]; FLOOR = float(sys.argv[5]) if len(sys.argv) > 5 else 50.0
os.makedirs(O, exist_ok=True)
R = pd.read_csv(mf, sep='\t', low_memory=False)
D = pd.read_csv(df, sep='\t')
G = pd.read_csv(gf, sep='\t')
R['exclude_code'] = ''; R['exclude_reason'] = ''
def mark(mask, code, reason):
    m = mask & (R.exclude_code == '')
    R.loc[m, 'exclude_code'] = code; R.loc[m, 'exclude_reason'] = reason
    return int(m.sum())

ov_lib = set(D[D.decision == 'override_lib'].id); ov_ena = set(D[D.decision == 'override_noena'].id)
for r in ov_ena: R.loc[R.run == r, 'BioProject'] = 'SRA_' + r          # NCBI-only run: its own study
mark(R.study_accession.isna() & ~R.run.isin(ov_ena), 'NOENA', 'run absent from ENA; library and sample metadata cannot be verified')
mark(R.library_strategy.notna() & (R.library_strategy != 'WGS') & ~R.BioProject.isin(ov_lib), 'LIB', 'library strategy ' + R.library_strategy.fillna('') + ' (not WGS shotgun)')
for _, d in D[D.decision == 'exclude'].iterrows():
    if d.level == 'run': m = R.run == d.id
    elif d.level == 'sample': m = R.sample_accession == d.id
    elif d.level == 'bioproject': m = R.BioProject == d.id
    elif d.level == 'run_list':
        bp, what = d.id.split(':')
        iso = R.isolation_source.fillna('').str.lower()
        m = (R.BioProject == bp) & (iso.str.contains('farm soil') if what == 'soil' else iso.str.contains('human gut'))
    elif d.level == 'bioproject_list':
        bp, what = d.id.split(':')
        m = (R.BioProject == bp) & R.country.fillna('').str.contains('Davies Reef')
    else: raise SystemExit('unknown level ' + d.level)
    n_any = int(m.sum())
    if n_any == 0: raise SystemExit(f'decision matched no run: {d.level} {d.id}')
    if (R.loc[m, 'host_system'] != d.host_system).any(): raise SystemExit(f'host mismatch for {d.id}')
    mark(m, d.code, d.reason)
mark(R.GE.isna(), 'NOGE', 'MicrobeCensus gave no genome-equivalent estimate (re-run also failed)')
pend = D[D.decision == 'pending']
R['pending_check'] = R.BioProject.isin(pend[pend.level == 'bioproject'].id)

# study groups
sg = {bp: g.study_group for _, g in G.iterrows() for bp in g.BioProjects.split(',')}
R['study'] = R.BioProject.map(sg).fillna(R.BioProject)

# BioSample pooling (included runs only)
R['sample_id'] = R.sample_accession.fillna(R.run)
inc = R[R.exclude_code == '']
S = inc.groupby('sample_id').agg(host_system=('host_system', 'first'), n_host=('host_system', 'nunique'), BioProject=('BioProject', 'first'),
                                 study=('study', 'first'), runs=('run', lambda s: ','.join(sorted(s))), n_runs=('run', 'size'),
                                 sets=('set', lambda s: ','.join(sorted(set(s)))), reads_processed=('reads_processed', 'sum'),
                                 GE=('GE', lambda s: s.sum(min_count=1) if s.notna().all() else np.nan), Gb=('Gb', 'sum'),
                                 reads_mapped=('reads_mapped', 'sum'), lat=('lat', 'first'), lon=('lon', 'first'), country=('country', 'first'),
                                 platform=('instrument_platform', 'first'), pending_check=('pending_check', 'max')).reset_index()
assert (S.n_host == 1).all(), 'a BioSample spans two host systems'
S['GE_per_M'] = S.GE / S.reads_processed * 1e6
S['set_main'] = S.GE > 0                 # main analysis: every valid sample (one method for all systems)
S['set_ge20'] = S.GE >= 20               # sensitivity analyses: minimum genome equivalents
S['set_ge50'] = S.GE >= 50
S['set_ge200'] = S.GE >= 200
S['set_full'] = S.GE > 0
R['in_main'] = R.sample_id.map(S.set_index('sample_id').set_main).fillna(False).astype(bool)
assert not ((R.exclude_code == '') & ~R.in_main).any()
R.to_csv(f'{O}/run_qc.tsv', sep='\t', index=False)
S.to_csv(f'{O}/samples.tsv', sep='\t', index=False)

H = ['Chicken', 'Coral', 'Rhizosphere', 'Sponge']
L = []
L.append('main set = all valid samples; GE >= 20 / 50 / 200 are sensitivity sets')
L.append('runs by exclusion code and host:\n' + R.groupby(['exclude_code', 'host_system']).size().unstack(fill_value=0).reindex(columns=H, fill_value=0).to_string())
for col in ['set_main', 'set_ge20', 'set_ge50', 'set_ge200']:
    x = S[S[col]]
    L.append(f'{col}: samples {len(x)} ' + str(x.host_system.value_counts().reindex(H).to_dict()) + ' | studies ' + str(x.groupby('host_system').study.nunique().reindex(H).to_dict()) + ' | BioProjects ' + str(x.groupby('host_system').BioProject.nunique().reindex(H).to_dict()))
L.append('pooled BioSamples (more than one run): ' + str(int((S.n_runs > 1).sum())) + ' samples from ' + str(int(S.loc[S.n_runs > 1, 'n_runs'].sum())) + ' runs')
L.append('samples still pending a metadata check (in main set): ' + str(S[S.set_main & S.pending_check].groupby(['host_system', 'BioProject']).size().to_dict()))
open(f'{O}/qc_summary.txt', 'w').write('\n\n'.join(L) + '\n'); print('\n\n'.join(L))
