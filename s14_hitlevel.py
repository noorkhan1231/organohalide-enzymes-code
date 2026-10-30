"""Step 14. Gene- and genome-level results (independent of the metagenome samples) on the 9,396 curated genes.
usage: python3 s14_hitlevel.py <inputs/hitlevel> <inputs/bundle> <out_dir>
Inputs (all on the curated set; computed earlier on the server and checked here):
  All_families_taxonomy_curated9396.tsv  hit -> genome, host system, family, GTDB lineage
  identity_to_reference_per_hit.csv      DIAMOND identity of each hit to its closest characterised HaloCycDB reference
  Domain_architecture_per_hit.xlsx       Pfam-A architecture per hit (hmmscan --cut_ga, i-Evalue <= 1e-5) + subfamily label
  Domain_vs_host_all_hits.xlsx, Domains_marine_terrestrial_and_pockets.xlsx   domain tests (host; marine vs terrestrial;
                                         within GTDB order, CMH) - copied, not recomputed
  deeptmhmm_per_hit_annotated.tsv, deeptmhmm_host_tests.tsv   DeepTMHMM topology per hit and host tests within order
  mob_hits_context.tsv, mob_background.tsv   geNomad context of Coral/Sponge RDase, HyDase, VanHase hits and background
  bundle/Denominators/genomes_<host>.txt     genomes searched per host system
Outputs: genes_per_family_host.csv, prevalence_copy_number.csv, identity_summary.csv, rdase_module.csv, vanhase_module.csv,
         deeptmhmm_summary.csv, mobility_plasmid_contig.csv, hitlevel_summary.json
"""
import sys, os, json, numpy as np, pandas as pd
from scipy.stats import fisher_exact, binomtest
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stats_utils import HOSTS, FAM, bh
H, B, O = sys.argv[1:4]; os.makedirs(O, exist_ok=True)
t = pd.read_csv(f'{H}/All_families_taxonomy_curated9396.tsv', sep='\t')
assert len(t) == 9396 and t.hit_id.is_unique, 'curated set must have 9,396 unique hits'
den = {h: sum(1 for _ in open(f'{B}/Denominators/genomes_{h}.txt')) for h in HOSTS}
S = dict(curated_genes=len(t), carrier_genomes=int(t.genome_id.nunique()), genomes_searched=sum(den.values()), searched_by_host=den,
         carriers_by_host=t.groupby('environment').genome_id.nunique().reindex(HOSTS).to_dict())
# genes per family x host and prevalence / copy number per carrier genome
g = t.groupby(['gene', 'environment']).agg(genes=('hit_id', 'size'), carriers=('genome_id', 'nunique')).reset_index()
g['mean_copies'] = g.genes / g.carriers
g['max_copies'] = [t[(t.gene == r.gene) & (t.environment == r.environment)].genome_id.value_counts().max() for r in g.itertuples()]
g['genomes_searched'] = g.environment.map(den); g['prevalence_pct_searched'] = 100 * g.carriers / g.genomes_searched
g['carriers_with_any_hit'] = g.environment.map(S['carriers_by_host'])
g['prevalence_pct_carriers'] = 100 * g.carriers / g.carriers_with_any_hit
g.rename(columns={'gene': 'family', 'environment': 'host_system'}).to_csv(f'{O}/prevalence_copy_number.csv', index=False)
pd.crosstab(t.gene, t.environment).reindex(index=FAM, columns=HOSTS, fill_value=0).to_csv(f'{O}/genes_per_family_host.csv')
# identity to closest characterised reference
i = pd.read_csv(f'{H}/identity_to_reference_per_hit.csv'); assert i.hit_id.isin(t.hit_id).all()
iq = i.groupby('family').identity.agg(n='size', median='median', q25=lambda v: v.quantile(.25), q75=lambda v: v.quantile(.75),
                                       pct_below_50=lambda v: 100 * (v < 50).mean()).reindex(FAM)
iq['hits'] = t.gene.value_counts().reindex(FAM); iq['no_alignment'] = iq.hits - iq.n
pooled = i[~i.family.isin(['GSTDase', 'MetDase'])]
iq.loc['pooled_excl_GSTDase_MetDase'] = [len(pooled), pooled.identity.median(), pooled.identity.quantile(.25), pooled.identity.quantile(.75),
                                         100 * (pooled.identity < 50).mean(), (~t.gene.isin(['GSTDase', 'MetDase'])).sum(), np.nan]
iq['no_alignment'] = iq.hits - iq.n
iq.to_csv(f'{O}/identity_summary.csv')
S['identity_pooled_excl_GSTDase_MetDase'] = dict(aligned=int(len(pooled)), median=float(pooled.identity.median()), pct_below_50=round(100 * (pooled.identity < 50).mean(), 1))
# domain architecture: RDase motif group electron-transfer module, VanHase extended N-terminal module
a = pd.read_excel(f'{H}/Domain_architecture_per_hit.xlsx'); assert set(a.hit_id) == set(t.hit_id)
r = a[(a.gene == 'RDase') & (a.subfamily_call == 'motifgroup_CxxCxxCxxxCP')].copy()
r['module'] = r.arch.fillna('').str.contains(r'(?:^|\|)(?:Fer2|DmmA-like_N)(?:\||$)', regex=True)
rt = r.groupby('environment').module.agg(['size', 'sum']).reindex(HOSTS).fillna(0).astype(int).rename(columns={'size': 'genes', 'sum': 'with_module'})
mar = rt.loc[['Coral', 'Sponge']].sum(); ter = rt.loc[['Chicken', 'Rhizosphere']].sum()
OR, P = fisher_exact([[mar.with_module, mar.genes - mar.with_module], [ter.with_module, ter.genes - ter.with_module]])
rt.to_csv(f'{O}/rdase_module.csv'); S['rdase_motifgroup_module'] = dict(marine=f'{mar.with_module}/{mar.genes}', terrestrial=f'{ter.with_module}/{ter.genes}', fisher_P=P)
v = a[a.gene == 'VanHase'].copy()
v['module'] = v.arch.fillna('').str.contains(r'FG-GAP|UnbV_ASPIC', regex=True)
vt = v.groupby(['environment']).module.agg(['size', 'sum']).reindex(HOSTS).fillna(0).astype(int).rename(columns={'size': 'genes', 'sum': 'with_module'})
vt.to_csv(f'{O}/vanhase_module.csv'); S['vanhase_FG-GAP_UnbV_module'] = vt.with_module.to_dict()
# DeepTMHMM: signal peptide per family x host (tests within GTDB order are copied from deeptmhmm_host_tests.tsv)
d = pd.read_csv(f'{H}/deeptmhmm_per_hit_annotated.tsv', sep='\t'); assert set(d.hit_id) == set(t.hit_id)
ds = d.groupby(['gene', 'environment']).agg(genes=('SP', 'size'), signal_peptide=('SP', 'sum'), tm_helix=('TMp', 'sum')).reset_index()
ds['SP_pct'] = 100 * ds.signal_peptide / ds.genes; ds.to_csv(f'{O}/deeptmhmm_summary.csv', index=False)
S['signal_peptide_total'] = int(d.SP.sum())
# geNomad: plasmid-classified contigs, expected = background rate of all genes in the same contig-length class
m = pd.read_csv(f'{H}/mob_hits_context.tsv', sep='\t'); bg = pd.read_csv(f'{H}/mob_background.tsv', sep='\t')
bw = bg.pivot_table(index=['environment', 'len_bin'], columns='metric', values='count')
rows = []
for (env, fam), x in m.groupby(['environment', 'family']):
    exp_p = sum(bw.loc[(env, lb), 'plasmid'] / bw.loc[(env, lb), 'genes'] * n for lb, n in x.len_bin.value_counts().items()) / len(x)
    k = int(x.plasmid_contig.sum()); p = binomtest(k, len(x), exp_p).pvalue
    exp_m = sum(bw.loc[(env, lb), 'mob10kb'] / bw.loc[(env, lb), 'genes'] * n for lb, n in x.len_bin.value_counts().items()) / len(x)
    km = int((x.n_mob_genes_10kb > 0).sum()); pm = binomtest(km, len(x), exp_m).pvalue
    rows.append(dict(host_system=env, family=fam, genes=len(x), plasmid_contig=k, plasmid_pct=100 * k / len(x), expected_pct=100 * exp_p, P=p,
                     mobility_gene_10kb=km, mobility_pct=100 * km / len(x), mobility_expected_pct=100 * exp_m, P_mobility=pm, in_provirus=int(x.in_provirus.sum())))
mb = pd.DataFrame(rows); fdr = bh(np.r_[mb.P.values, mb.P_mobility.values]); mb['FDR'] = fdr[:len(mb)]; mb['FDR_mobility'] = fdr[len(mb):]
mb.to_csv(f'{O}/mobility_plasmid_contig.csv', index=False)
for f in ['Domain_vs_host_all_hits.xlsx', 'Domains_marine_terrestrial_and_pockets.xlsx', 'deeptmhmm_host_tests.tsv']:
    os.system(f'cp "{H}/{f}" "{O}/"')
json.dump(S, open(f'{O}/hitlevel_summary.json', 'w'), indent=1, default=float)
pd.set_option('display.width', 220)
print(json.dumps(S, indent=1, default=float)); print(g.round(2).to_string(index=False)); print(iq.round(1).to_string()); print(rt, '\n', vt)
print(mb.round(4).to_string(index=False))
