"""Step 7. Sample-level abundance tables: RPKG = reads per kilobase of target per genome equivalent, all targets.
usage: python3 s04_build.py <inputs_dir> <s03 dir> <out_dir>
Counts of runs belonging to one BioSample are summed, genome equivalents are summed (samples.tsv), so each BioSample
is one observation. RPK_t = NumReads_t / (Length_t / 1000); RPKG(sample, set) = sum_t RPK_t / GE(sample).
Writes, for every analysis set (full, ge20, ge50, main, ge200): family_rpkg_<set>.csv, subfamily_rpkg_long_<set>.csv, and
genus_subfamily_rpkg_main.csv; plus native-host sensitivity tables (targets from the sample's own host system only).
"""
import sys, os, glob, pandas as pd, numpy as np
I, Q, O = sys.argv[1:4]; os.makedirs(O, exist_ok=True)
FAM = ['RDase', 'HyDase', 'GSTDase', 'OxDase', 'DehyDase', 'MetDase', 'FlaHase', 'VanHase', 'NHFeHase', 'SAMHase']
rq = pd.read_csv(f'{Q}/run_qc.tsv', sep='\t', low_memory=False)
S = pd.read_csv(f'{Q}/samples.tsv', sep='\t')
meta = pd.read_csv(f'{I}/original/target_meta.tsv', sep='\t')
files = [f'{I}/original/per_target_counts.tsv.gz'] + sorted(glob.glob(f'{I}/job*/extra_counts.tsv.gz'))
c = pd.concat([pd.read_csv(f, sep='\t', dtype={'run': str}) for f in files], ignore_index=True)
c = c[c.dataset.isin(['host', 'extra'])]
r2s = rq[rq.exclude_code.isin(['', 'LOWGE']) | rq.exclude_code.isna()].set_index('run').sample_id
c = c[c.run.isin(r2s.index)].copy(); c['sample_id'] = c.run.map(r2s)
dup = c.duplicated(['run', 'family', 'target']).sum(); assert dup == 0, f'{dup} duplicated count rows'
m = c.merge(meta, on=['target', 'family'], how='left'); assert m.origin_host.notna().all(), 'targets without metadata'
m = m.groupby(['sample_id', 'family', 'target', 'Length', 'origin_host', 'genus', 'order', 'subfamily'], dropna=False).NumReads.sum().reset_index()
m['RPK'] = m.NumReads / (m.Length / 1000.0)
Si = S.set_index('sample_id')
m['host_system'] = m.sample_id.map(Si.host_system); m['GE'] = m.sample_id.map(Si.GE)
m = m[m.GE > 0]; m['RPKG'] = m.RPK / m.GE; m['native'] = m.origin_host == m.host_system
print('count rows', len(m), '| samples with any read', m.sample_id.nunique(), '| reads', int(m.NumReads.sum()))

def wide(df, col, ids):
    w = df.groupby(['sample_id', 'family'])[col].sum().unstack(fill_value=0).reindex(columns=FAM, fill_value=0)
    base = Si.loc[ids, ['host_system', 'study', 'BioProject', 'GE', 'reads_processed', 'platform']]
    return base.join(w, how='left').fillna({f: 0 for f in FAM})

for s in ['full', 'ge20', 'ge50', 'main', 'ge200']:
    ids = S.loc[S[f'set_{s}'], 'sample_id']
    wide(m, 'RPKG', ids).to_csv(f'{O}/family_rpkg_{s}.csv')
    if s == 'main':
        wide(m[m.native], 'RPKG', ids).to_csv(f'{O}/family_rpkg_native_main.csv')
        wide(m, 'NumReads', ids).to_csv(f'{O}/family_reads_main.csv')
        g = m[m.sample_id.isin(ids)].groupby(['family', 'subfamily', 'genus', 'origin_host']).agg(RPKG=('RPKG', 'sum'), reads=('NumReads', 'sum'), n_samples=('sample_id', 'nunique')).reset_index()
        g.to_csv(f'{O}/genus_subfamily_rpkg_main.csv', index=False)
    sub = m[m.sample_id.isin(ids)].groupby(['sample_id', 'family', 'subfamily']).agg(RPKG=('RPKG', 'sum'), reads=('NumReads', 'sum')).reset_index()
    sub.to_csv(f'{O}/subfamily_rpkg_long_{s}.csv', index=False)
    w = pd.read_csv(f'{O}/family_rpkg_{s}.csv')
    print(s, 'samples', len(w), w.host_system.value_counts().to_dict(), '| all-zero', int((w[FAM].sum(axis=1) == 0).sum()))
ns = m[m.sample_id.isin(S.loc[S.set_main, 'sample_id'])].groupby(['host_system', 'family']).apply(lambda x: x.loc[x.native, 'NumReads'].sum() / x.NumReads.sum(), include_groups=False).unstack()
ns.to_csv(f'{O}/native_read_share_main.csv'); print('native read share (main):'); print(ns.round(3).to_string())
