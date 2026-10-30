"""Step 11. Host-associated vs surrounding (ambient) metagenomes, at study level, with a system-matched null-gene control.
usage: python3 s11_ambient.py <inputs_dir> <s03_dir> <s04_dir> <out_dir>

Host side = the final dataset (main set; runs that passed QC, s03 run_qc.tsv). Ambient side = seawater (Coral, Sponge), bulk
soil (Rhizosphere) and poultry-house / farm environment (Chicken) metagenomes, processed with the same fastp + Salmon settings
(inputs/bundle/Ambient lists; inputs/original per_target_counts / genome_equivalents / quant_runs, dataset 'ambient' or 'host').

Null-gene control: Salmon on the first 10 million read pairs against random non-target genes drawn from the genome
collection of each host system (inputs/original/null_summary.tsv = earlier host runs + all ambient runs;
inputs/null_extra/null_summary_extra.tsv = runs added for the remaining host studies, 2-3 runs per study).
null RPKG(run, system) = null RPK on the subsample x (reads processed at full depth / reads processed in subsample) / GE(run).
Null-normalised abundance(run) = RPKG(family, run) / null RPKG(run, matching system); it removes the advantage a host sample
has because the reference genomes were assembled from host-associated samples.

Tests (per host system, per family): two-sided Mann-Whitney U on study means (host: study = BioProject or merged study
group used in all other analyses; ambient: BioProject), only when both sides have >= 3 studies; BH FDR within
scope x measure. Measures: RPKG (all main-set host samples, sample level from s04) and null_normalised (runs with a null value).
Call (conservative): a family is host-enriched (or surroundings-enriched) only when RPKG and null-normalised abundance both
differ in the same direction at FDR < 0.05; opposite directions = 'unresolved'. Matched-only comparisons (ambient runs from
the host BioProjects) are not tested: only Sponge reaches 3 matched studies.
"""
import sys, os, csv, json, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stats_utils import HOSTS, FAM, bh
from scipy.stats import mannwhitneyu
I, Q, B, out = sys.argv[1:5]; os.makedirs(out, exist_ok=True)
A = f'{I}/bundle/Ambient'

# ---- ambient run list (unchanged from the earlier ambient analysis)
EXT = {'seawater': ['Coral', 'Sponge'], 'soil': ['Rhizosphere'], 'chicken': ['Chicken']}
amb = {}
for r in csv.DictReader(open(f'{A}/ambient_matched_own_projects.tsv'), delimiter='\t'): amb[r['run_accession']] = ([r['our_system']], 'matched', r['BioProject'])
for f, k in [('seawater_external.tsv', 'seawater'), ('bulk_soil_external.tsv', 'soil'), ('chicken_environment_external.tsv', 'chicken')]:
    for r in csv.DictReader(open(f'{A}/{f}'), delimiter='\t'): amb[r['Run']] = (EXT[k], 'external', r['BioProject'])
for r in csv.DictReader(open(f'{A}/topup_candidates.tsv'), delimiter='\t'): amb[r['Run']] = (EXT[r['category']], 'external', r['BioProject'])
for x in ('SRR30894617', 'SRR30894618', 'SRR30894619'): amb.pop(x, None)          # withdrawn at NCBI
ro = pd.read_csv(f'{Q}/runs_original_for_ambient.tsv', sep='\t')
for r in ro[ro.ambient_for.notna()].itertuples(): amb[r.run] = ([r.ambient_for], 'matched', r.BioProject)   # farm soil of a chicken study

# ---- host runs in the final dataset
rq = pd.read_csv(f'{Q}/run_qc.tsv', sep='\t', low_memory=False)
rq = rq[rq.exclude_code.isna() | rq.exclude_code.isin(['', 'LOWGE'])]
S_ = pd.read_csv(f'{Q}/samples.tsv', sep='\t'); main_ids = set(S_.loc[S_.set_main, 'sample_id'])
rq = rq[rq.sample_id.isin(main_ids)]
clash = [r for r in amb if r in set(rq.run)]
for r in clash: amb.pop(r)                                    # a run cannot be host and ambient at the same time
print('ambient runs that are host runs in the final dataset (removed from ambient):', len(clash))

# ---- run-level counts, full-depth reads processed, GE
cnt = [pd.read_csv(f'{I}/original/per_target_counts.tsv.gz', sep='\t', dtype={'run': str})]
qr = [pd.read_csv(f'{I}/original/quant_runs.tsv', sep='\t', dtype={'run': str})]
for j in ['job6', 'job6b', 'job14', 'job16']:
    cnt.append(pd.read_csv(f'{I}/{j}/extra_counts.tsv.gz', sep='\t', dtype={'run': str}))
    qr.append(pd.read_csv(f'{I}/{j}/extra_quant_runs.tsv', sep='\t', dtype={'run': str}))
c = pd.concat(cnt, ignore_index=True); q = pd.concat(qr, ignore_index=True)
c['ds'] = np.where(c.dataset == 'ambient', 'ambient', 'host'); q['ds'] = np.where(q.dataset == 'ambient', 'ambient', 'host')
c = c[c.dataset != 'metatx']; q = q[q.dataset != 'metatx']
c['RPK'] = c.NumReads / (c.Length / 1000)
fam = c.groupby(['ds', 'run', 'family']).RPK.sum().unstack(fill_value=0).reindex(columns=FAM, fill_value=0)
nproc = q.groupby(['ds', 'run']).num_processed.max()
ge0 = pd.read_csv(f'{I}/original/genome_equivalents.tsv', sep='\t').drop_duplicates(['dataset', 'run']).set_index(['dataset', 'run']).genome_equivalents

rows = []
for r in rq.itertuples(): rows.append(dict(ds='host', run=r.run, system=r.host_system, group='host', unit=r.study, sample_id=r.sample_id, GE=r.GE))
for run, (systems, kind, bp) in amb.items():
    ds = 'ambient' if ('ambient', run) in ge0.index else 'host'
    for s in systems: rows.append(dict(ds=ds, run=run, system=s, group='ambient_' + kind, unit=bp, sample_id=run, GE=ge0.get((ds, run), np.nan)))
R = pd.DataFrame(rows)
R['nproc_full'] = [nproc.get((d, r), np.nan) for d, r in zip(R.ds, R.run)]
R = R.join(fam, on=['ds', 'run'])
miss_q = R[R.nproc_full.isna()]
print('runs without full-depth quant (dropped):', len(miss_q), miss_q.groupby('group').size().to_dict())
R = R[R.nproc_full.notna() & (R.GE > 0)].fillna({f: 0 for f in FAM})
for f in FAM: R[f] = R[f] / R.GE

# ---- null genes
n0 = pd.read_csv(f'{I}/original/null_summary.tsv', sep='\t'); n1 = pd.read_csv(f'{I}/null_extra/null_summary_extra.tsv', sep='\t')
n0['src'] = 'earlier'; n1['src'] = 'null_extra'
nl = pd.concat([n0, n1], ignore_index=True); nl['ds'] = np.where(nl.group == 'ambient', 'ambient', 'host')
assert not nl.duplicated(['ds', 'run', 'gene_system']).any(), 'duplicated null rows'
nrpk = nl.set_index(['ds', 'run', 'gene_system']).rpk; nsub = nl.groupby(['ds', 'run']).num_processed_subsample.max()
nsrc = nl.groupby(['ds', 'run']).src.first()
R['null_rpk_sub'] = [nrpk.get((d, r, s), np.nan) for d, r, s in zip(R.ds, R.run, R.system)]
R['nproc_sub'] = [nsub.get((d, r), np.nan) for d, r in zip(R.ds, R.run)]
R['null_source'] = [nsrc.get((d, r), '') for d, r in zip(R.ds, R.run)]
R['null_RPKG'] = R.null_rpk_sub * (R.nproc_full / R.nproc_sub) / R.GE
R.to_csv(f'{out}/ambient_runs_rpkg.csv', index=False)

# host RPKG at sample level (all main-set samples, s04), ambient at run level (each ambient run = one BioSample)
fs = pd.read_csv(f'{B}/family_rpkg_main.csv')
H = fs[['sample_id', 'host_system', 'study'] + FAM].rename(columns={'host_system': 'system', 'study': 'unit'}).assign(group='host')
Am = R[R.group != 'host'][['sample_id', 'system', 'unit', 'group'] + FAM]
SAMP = pd.concat([H, Am], ignore_index=True)

def proj_test(x, col):
    h = x[x.group == 'host'].groupby('unit')[col].mean(); a = x[x.group != 'host'].groupby('unit')[col].mean()
    d = dict(n_host_studies=len(h), n_ambient_studies=len(a), n_host=int((x.group == 'host').sum()), n_ambient=int((x.group != 'host').sum()))
    if len(h) < 3 or len(a) < 3: return d
    p = mannwhitneyu(h, a, alternative='two-sided').pvalue
    return dict(d, mean_host=h.mean(), mean_ambient=a.mean(), median_host=h.median(), median_ambient=a.median(),
                ratio=h.mean() / a.mean() if a.mean() > 0 else (np.inf if h.mean() > 0 else np.nan),
                host_studies_detected=int((h > 0).sum()), ambient_studies_detected=int((a > 0).sum()), P=p)
T = []
for sysn in HOSTS:
    xs = SAMP[SAMP.system == sysn]; xr = R[(R.system == sysn) & (R.null_RPKG > 0)].copy()
    for scope in ['all_ambient']:
        ks = xs.group.notna(); kr = xr.group.notna()
        T.append(dict(system=sysn, scope=scope, set='NULL', measure='RPKG', **proj_test(xr[kr], 'null_RPKG')))
        for f in FAM:
            T.append(dict(system=sysn, scope=scope, set=f, measure='RPKG', **proj_test(xs[ks], f)))
            x2 = xr[kr].copy(); x2['norm'] = x2[f] / x2.null_RPKG
            T.append(dict(system=sysn, scope=scope, set=f, measure='null_normalised', **proj_test(x2, 'norm')))
T = pd.DataFrame(T); T['FDR'] = np.nan
for k, g in T.groupby(['scope', 'measure']):
    T.loc[g.index, 'FDR'] = bh(g.P.values)
T.to_csv(f'{out}/host_vs_ambient_family.csv', index=False)

cov = R.groupby(['system', 'group']).agg(runs=('run', 'nunique'), studies=('unit', 'nunique'),
                                         runs_with_null=('null_RPKG', lambda v: int((v > 0).sum())))
cov['studies_with_null'] = R[R.null_RPKG > 0].groupby(['system', 'group']).unit.nunique()
cov.to_csv(f'{out}/ambient_coverage.csv')
pd.set_option('display.width', 250)
print(cov.to_string())
print(T[T.scope == 'all_ambient'].round(4).to_string(index=False))
a = T[T.measure == 'RPKG'].set_index(['system', 'set']); b = T[T.measure == 'null_normalised'].set_index(['system', 'set'])
C = a[['n_host_studies', 'n_ambient_studies', 'host_studies_detected', 'ambient_studies_detected', 'ratio', 'P', 'FDR']].join(
    b[['ratio', 'P', 'FDR']], rsuffix='_null_normalised').reset_index()
C = C[C.set != 'NULL']
def call(r):
    if not (r.FDR < 0.05 and r.FDR_null_normalised < 0.05): return 'no consistent difference'
    if r.ratio > 1 and r.ratio_null_normalised > 1: return 'host-enriched'
    if r.ratio < 1 and r.ratio_null_normalised < 1: return 'surroundings-enriched'
    return 'unresolved'
C['call'] = C.apply(call, axis=1); C.to_csv(f'{out}/host_vs_ambient_calls.csv', index=False)
print(C.round(4).to_string(index=False))
json.dump({s: C[C.system == s].groupby('call').set.apply(list).to_dict() for s in HOSTS}, open(f'{out}/ambient_summary.json', 'w'), indent=1)
null = T[T.set == 'NULL'][['system', 'n_host_studies', 'n_ambient_studies', 'mean_host', 'mean_ambient', 'ratio', 'P', 'FDR']]
null.to_csv(f'{out}/null_gene_recruitment.csv', index=False); print(null.round(4).to_string(index=False))
