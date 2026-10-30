"""Step 4. One metadata + per-run QC table for every metagenome run (original + jobs 6, 6b, 14, 16).
usage: python3 s02_metadata.py <inputs_dir> <out_dir>
Columns: run, set, host_system, BioProject, ENA fields (sample, titles, host, isolation source, lat/lon, country,
library, platform, bases), Salmon reads processed, genome equivalents (MicrobeCensus), fastp before/after reads and
bases, Q30, mean read length, targets detected, reads mapped, mapped reads per million, genome equivalents per million reads.
"""
import sys, os, glob, json, pandas as pd, numpy as np
I, O = sys.argv[1:3]; os.makedirs(O, exist_ok=True)
SETS = {'original': ('original/quant_runs.tsv', 'original/genome_equivalents.tsv', 'original/per_target_counts.tsv.gz', None),
        'job6': ('job6/extra_quant_runs.tsv', 'job6/extra_ge.tsv', 'job6/extra_counts.tsv.gz', 'job6/fastp_json'),
        'job6b': ('job6b/extra_quant_runs.tsv', 'job6b/extra_ge.tsv', 'job6b/extra_counts.tsv.gz', 'job6b/fastp_json'),
        'job14': ('job14/extra_quant_runs.tsv', 'job14/extra_ge.tsv', 'job14/extra_counts.tsv.gz', 'job14/fastp_json'),
        'job16': ('job16/extra_quant_runs.tsv', 'job16/extra_ge.tsv', 'job16/extra_counts.tsv.gz', 'job16/fastp_json')}
GE_FIX = f'{I}/ge_fix.tsv'   # run, genome_equivalents from the MicrobeCensus re-runs (optional)
rows = []
for s, (qf, gf, cf, fj) in SETS.items():
    if not os.path.exists(f'{I}/{qf}'): print('not available yet:', s); continue
    q = pd.read_csv(f'{I}/{qf}', sep='\t'); q = q[q.dataset.isin(['host', 'extra'])]
    g = pd.read_csv(f'{I}/{gf}', sep='\t'); g = g[g.dataset.isin(['host', 'extra'])]
    g['genome_equivalents'] = pd.to_numeric(g.genome_equivalents, errors='coerce')
    c = pd.read_csv(f'{I}/{cf}', sep='\t'); c = c[c.dataset.isin(['host', 'extra'])]
    r = q.groupby('run').agg(host_system=('host_dir', 'first'), reads_processed=('num_processed', 'max')).reset_index()
    r['set'] = s
    r['GE'] = r.run.map(g.dropna(subset=['genome_equivalents']).drop_duplicates('run').set_index('run').genome_equivalents)
    cc = c.groupby('run').agg(targets_detected=('target', 'nunique'), reads_mapped=('NumReads', 'sum'))
    r = r.join(cc, on='run'); r[['targets_detected', 'reads_mapped']] = r[['targets_detected', 'reads_mapped']].fillna(0)
    if fj:
        fq = []
        for f in glob.glob(f'{I}/{fj}/*.json'):
            try: d = json.load(open(f))['summary']
            except Exception: continue
            b, a = d['before_filtering'], d['after_filtering']
            fq.append(dict(run=os.path.basename(f)[:-5], fastp_reads_in=b['total_reads'], fastp_bases_in=b['total_bases'],
                           fastp_reads_out=a['total_reads'], fastp_bases_out=a['total_bases'], q30_after=a['q30_rate'],
                           mean_len_after=a['read1_mean_length'], gc_after=a['gc_content']))
        if fq: r = r.merge(pd.DataFrame(fq).drop_duplicates('run'), on='run', how='left')
    rows.append(r)
R = pd.concat(rows, ignore_index=True)
if os.path.exists(GE_FIX):
    fx = pd.read_csv(GE_FIX, sep='\t').set_index('run').genome_equivalents
    m = R.GE.isna() & R.run.isin(fx.index); R.loc[m, 'GE'] = R.loc[m, 'run'].map(fx); print('GE filled from re-runs:', int(m.sum()))
R['mapped_per_M'] = R.reads_mapped / R.reads_processed * 1e6
R['GE_per_M'] = R.GE / R.reads_processed * 1e6

ena = pd.read_csv(f'{I}/metadata/ena_metadata_all_runs.tsv', sep='\t', dtype=str)
for extra in glob.glob(f'{I}/metadata/ena_job6b.tsv'):
    ena = pd.concat([ena, pd.read_csv(extra, sep='\t', dtype=str)], ignore_index=True)
ena = ena.drop_duplicates('run_accession').rename(columns={'run_accession': 'run'})
R = R.merge(ena, on='run', how='left')
R['BioProject'] = R.study_accession
# original runs: BioProject from the curated run table (identical to ENA study accession; asserted)
orig = pd.read_csv(f'{I}/original/runs.tsv', sep='\t').set_index('run')
o = R.set == 'original'
mism = R[o & (R.run.map(orig.BioProject) != R.BioProject)]
print('original runs whose ENA study differs from curated BioProject:', len(mism), mism[['run', 'BioProject']].assign(cur=mism.run.map(orig.BioProject)).values.tolist()[:10])
R.loc[o, 'BioProject'] = R.loc[o, 'run'].map(orig.BioProject)
R['base_count'] = pd.to_numeric(R.base_count, errors='coerce'); R['Gb'] = R.base_count / 1e9
for k in ['lat', 'lon']: R[k] = pd.to_numeric(R[k], errors='coerce')
print('runs', len(R), R.groupby(['set', 'host_system']).size().to_dict())
print('missing ENA metadata:', int(R.study_accession.isna().sum()), '| missing GE:', int(R.GE.isna().sum()), R[R.GE.isna()].run.tolist())
R.to_csv(f'{O}/all_runs_metadata_qc.tsv', sep='\t', index=False)
