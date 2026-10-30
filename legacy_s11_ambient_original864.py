"""Step 11. Host-associated vs surrounding (ambient) metagenomes, at BioProject level, with a system-matched null-gene control.
usage: python3 05_ambient.py <counts.tsv.gz> <quant_runs.tsv> <genome_equivalents.tsv> <target_meta.tsv> <runs.tsv> <ambient_dir> <null_summary.tsv> <outdir>
Each family (and labelled subfamily) is expressed as RPKG; the null is RPKG of random non-target genes drawn from the same
genome collection as the host system (scaled from the 10 M read-pair subsample to full depth with Salmon num_processed).
Null-normalised abundance = RPKG(target set) / RPKG(null genes of that system): it removes the advantage a host sample
has because the reference genomes were assembled from host samples."""
import sys, os, csv, json, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stats_utils import *
from scipy.stats import mannwhitneyu
cf, qf, gef, mf, rf, AMB, nf, out = sys.argv[1:9]; os.makedirs(out, exist_ok=True)
EXT = {'seawater': ['Coral', 'Sponge'], 'soil': ['Rhizosphere'], 'chicken': ['Chicken']}
amb = {}
for r in csv.DictReader(open(f'{AMB}/ambient_matched_own_projects.tsv'), delimiter='\t'): amb[r['run_accession']] = ([r['our_system']], 'matched', r['BioProject'])
for f, k in [('seawater_external.tsv', 'seawater'), ('bulk_soil_external.tsv', 'soil'), ('chicken_environment_external.tsv', 'chicken')]:
    for r in csv.DictReader(open(f'{AMB}/{f}'), delimiter='\t'): amb[r['Run']] = (EXT[k], 'external', r['BioProject'])
for r in csv.DictReader(open(f'{AMB}/topup_candidates.tsv'), delimiter='\t'): amb[r['Run']] = (EXT[r['category']], 'external', r['BioProject'])
for x in ('SRR30894617', 'SRR30894618', 'SRR30894619'): amb.pop(x, None)          # withdrawn at NCBI
runs = pd.read_csv(rf, sep='\t')
extra = runs[runs.get('ambient_for', pd.Series(index=runs.index, dtype=object)).notna()] if 'ambient_for' in runs.columns else runs.head(0)
for r in extra.itertuples(): amb[r.run] = ([r.ambient_for], 'matched', r.BioProject)   # e.g. farm soil of a chicken study
c = pd.read_csv(cf, sep='\t'); meta = pd.read_csv(mf, sep='\t')
c = c.merge(meta[['target', 'family', 'subfamily']], on=['target', 'family'], how='left')
c['RPK'] = c.NumReads / (c.Length / 1000)
ge = pd.read_csv(gef, sep='\t').drop_duplicates(['dataset', 'run']).set_index(['dataset', 'run']).genome_equivalents
q = pd.read_csv(qf, sep='\t').groupby(['dataset', 'run']).num_processed.max()
host_runs = runs[runs.include == 1].set_index('run')
rows = []
def add(sel_ds, run, systems, grp, bp):
    for sysn in systems: rows.append(dict(dataset=sel_ds, run=run, system=sysn, group=grp, BioProject=bp))
for r in host_runs.itertuples(): add('host', r.Index, [r.host_system], 'host', r.BioProject)
for run, (systems, kind, bp) in amb.items():
    ds = 'host' if run in runs.run.values else 'ambient'
    add(ds, run, systems, 'ambient_' + kind, bp)
S = pd.DataFrame(rows)
S['GE'] = [ge.get((d, r), np.nan) for d, r in zip(S.dataset, S.run)]
S['nproc_full'] = [q.get((d, r), np.nan) for d, r in zip(S.dataset, S.run)]
fam = c.groupby(['dataset', 'run', 'family']).RPK.sum().unstack(fill_value=0).reindex(columns=FAM, fill_value=0)
S = S.join(fam, on=['dataset', 'run']).fillna({f: 0 for f in FAM})
nl = pd.read_csv(nf, sep='\t'); nl['dataset'] = nl.group
nl = nl.pivot_table(index=['dataset', 'run'], columns='gene_system', values='rpk', aggfunc='sum')
nsub = pd.read_csv(nf, sep='\t').groupby(['group', 'run']).num_processed_subsample.max()
S['null_rpk_sub'] = [nl.loc[(d, r), s] if (d, r) in nl.index and s in nl.columns else np.nan for d, r, s in zip(S.dataset, S.run, S.system)]
S['nproc_sub'] = [nsub.get((d, r), np.nan) for d, r in zip(S.dataset, S.run)]
S['null_RPKG'] = S.null_rpk_sub * (S.nproc_full / S.nproc_sub) / S.GE
for f in FAM: S[f] = S[f] / S.GE
S = S[S.GE.notna()]
S.to_csv(f'{out}/ambient_samples_rpkg.csv', index=False)
def proj_test(x, col, S_):
    h = x[x.group == 'host'].groupby('BioProject')[col].mean(); a = x[x.group != 'host'].groupby('BioProject')[col].mean()
    if len(h) < 3 or len(a) < 3: return dict(n_host_proj=len(h), n_amb_proj=len(a))
    p = mannwhitneyu(h, a, alternative='two-sided').pvalue
    return dict(n_host_proj=len(h), n_amb_proj=len(a), mean_host=h.mean(), mean_ambient=a.mean(), ratio=(h.mean() + 1e-12) / (a.mean() + 1e-12), P=p)
T = []
for sysn in HOSTS:
    x = S[S.system == sysn].copy()
    for scope, xx in [('all_ambient', x), ('matched_only', x[(x.group == 'host') | (x.group == 'ambient_matched')])]:
        T.append(dict(system=sysn, scope=scope, set='NULL', measure='RPKG', **proj_test(xx, 'null_RPKG', sysn)))
        for f in FAM:
            T.append(dict(system=sysn, scope=scope, set=f, measure='RPKG', **proj_test(xx, f, sysn)))
            xx2 = xx[xx.null_RPKG > 0].copy(); xx2['norm'] = xx2[f] / xx2.null_RPKG
            T.append(dict(system=sysn, scope=scope, set=f, measure='null_normalised', **proj_test(xx2, 'norm', sysn)))
T = pd.DataFrame(T); T['FDR'] = np.nan
for k, g in T.groupby(['scope', 'measure']):
    T.loc[g.index, 'FDR'] = bh(g.P.values)
T.to_csv(f'{out}/host_vs_ambient_family.csv', index=False)
pd.set_option('display.width', 250)
print(S.groupby(['system', 'group']).agg(samples=('run', 'nunique'), projects=('BioProject', 'nunique')).to_string())
print(T[(T.scope == 'all_ambient')].round(4).to_string(index=False))
