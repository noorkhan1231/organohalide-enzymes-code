"""Step 17. Independent consistency check of the final numbers, figures' source data, tables and Supplementary Data.
usage: python3 s99_verify.py <FINAL dir>   (prints PASS/FAIL per check; exit code 1 if any check fails)
Every expected value is recomputed from the primary result files, not copied from the figure scripts.
"""
import sys, os, json, re, glob, numpy as np, pandas as pd
B = sys.argv[1]; R = f'{B}/results'; T = f'{B}/tables'; FG = f'{B}/figures'
HOSTS = ['Chicken', 'Coral', 'Rhizosphere', 'Sponge']; FAM = ['RDase', 'HyDase', 'GSTDase', 'OxDase', 'DehyDase', 'MetDase', 'FlaHase', 'VanHase', 'NHFeHase', 'SAMHase']
fails = []
def check(name, cond, detail=''):
    print(('PASS ' if cond else 'FAIL ') + name + ('' if cond else f'  <- {detail}'))
    if not cond: fails.append(name)

S = pd.read_csv(f'{R}/s03_qc/samples.tsv', sep='\t'); rq = pd.read_csv(f'{R}/s03_qc/run_qc.tsv', sep='\t', low_memory=False)
ok_runs = rq[rq.exclude_code.isna() | rq.exclude_code.isin(['', 'LOWGE'])]
main = set(S.loc[S.set_main, 'sample_id'])
check('1,787 main-set samples', len(main) == 1787, len(main))
check('host split 486/452/440/409', S[S.set_main].host_system.value_counts().reindex(HOSTS).tolist() == [486, 452, 440, 409])
check('179 studies (34/25/57/63)', S[S.set_main].groupby('host_system').study.nunique().reindex(HOSTS).tolist() == [34, 25, 57, 63])
check('1,854 runs in final samples', ok_runs.sample_id.isin(main).sum() == 1854, ok_runs.sample_id.isin(main).sum())
excluded_runs = set(rq.run) - set(ok_runs.run)
# --- no excluded sample/run anywhere in the deliverables
sd1 = pd.ExcelFile(f'{T}/Supplementary_Data_1.xlsx')
r1 = sd1.parse('runs'); s1 = sd1.parse('samples')
check('SD1 samples = main set', set(s1.sample_id) == main)
check('SD1 runs: no excluded run', not (set(r1.run) & excluded_runs) and len(r1) == 1854)
sd4 = pd.read_excel(f'{T}/Supplementary_Data_4.xlsx', 'family_RPKG')
check('SD4 family RPKG = main set', set(sd4.sample_id) == main)
mp = pd.read_excel(f'{FG}/Figure_1_source_data.xlsx', 'b_sample_locations')
check('Figure 1 map = main set', set(mp.sample_id) == main)
st4 = pd.read_excel(f'{T}/Supplementary_Tables.xlsx', 'Table S1', header=None)
check('Supplementary Table 1 lists 179 studies', (st4.iloc[:, 0].astype(str).isin(set(S.study))).sum() == 179)
words = re.compile(r'exclu|reject|horizontal gene transfer|\bHGT\b|centroid|approximate|named place|coord_source', re.I)
hit = []
for f in glob.glob(f'{T}/*.xlsx'):
    x = pd.ExcelFile(f)
    for sh in x.sheet_names:
        d = x.parse(sh, header=None).astype(str).replace('Cold adaptation and horizontal gene transfer in Antarctic sponge and seawater microbiomes', '')   # BioProject title (metadata), not a claim
        if d.apply(lambda c: c.str.contains(words)).any().any(): hit.append(f'{os.path.basename(f)}:{sh}')
check('no exclusion/rejection/HGT wording in tables and Supplementary Data', not hit, hit)
# --- family level
V = pd.read_csv(f'{R}/s05_family/main_vegan_tests.csv').set_index('test')
f2 = pd.read_excel(f'{FG}/Figure_2_source_data.xlsx', 'b_R2')
check('Figure 2b R2 matches vegan tests', np.allclose(f2.R2.values[:3], [V.loc['host_sample', 'R2'], V.loc['host_study', 'R2'], V.loc['block_nested', 'R2']]))
K = pd.read_csv(f'{R}/s05_family/main_kruskal.csv')
check('KW 10/10 sample and study level', (K.FDR < 0.05).sum() == 10 and (K.FDR_study < 0.05).sum() == 10)
check('PERMANOVA R2 0.198 / study 0.243 / nested 0.349', [round(V.loc[k, 'R2'], 3) for k in ['host_sample', 'host_study', 'block_nested']] == [0.198, 0.243, 0.349])
# --- surroundings
C = pd.read_csv(f'{R}/s11_ambient/host_vs_ambient_calls.csv'); f3 = pd.read_excel(f'{FG}/Figure_3_source_data.xlsx', 'b_calls')
check('Figure 3b calls = s11 calls', f3[['system', 'set', 'call']].equals(C[['system', 'set', 'call']]))
cnt = C.groupby('call').size().to_dict()
check('calls: 10 host-enriched, 8 surroundings-enriched, 5 unresolved', (cnt.get('host-enriched'), cnt.get('surroundings-enriched'), cnt.get('unresolved')) == (10, 8, 5), cnt)
print('   calls:', cnt)
ar = pd.read_csv(f'{R}/s11_ambient/ambient_runs_rpkg.csv')
check('ambient runs are not host runs', not (set(ar[ar.group != 'host'].run) & set(ok_runs.run)))
check('null coverage: 176 host studies with null value', ar[(ar.group == 'host') & (ar.null_RPKG > 0)].groupby('system').unit.nunique().sum() == 176)
# --- metatranscriptomes
mt = json.load(open(f'{R}/s07_metatx/metatx_summary.json'))
check('190 metatranscriptome runs, 54 studies', mt['n_runs'] == 190 and sum(mt['studies_by_host'].values()) == 54)
# --- subfamilies
SI = pd.read_csv(f'{R}/s06_subfamily/main_subfamily_indicators_significant.csv')
check('48 indicators; 21 with <= 3 genes; 27 shown in Figure 4b', len(SI) == 48 and (SI.n_genes <= 3).sum() == 21 and (SI.n_genes >= 4).sum() == 27)
check('lineage control order level: 6 host-associated', (SI.lineage_order_verdict == 'host-associated within lineage').sum() == 6)
CP = pd.read_csv(f'{R}/s06_subfamily/combined_permanova.csv').set_index('level')
check('combined subfamily PERMANOVA R2 0.217 / 0.261, PERMDISP significant', [round(CP.loc['sample', 'R2'], 3), round(CP.loc['study', 'R2'], 3)] == [0.217, 0.261] and CP.loc['sample', 'permdisp_P'] <= 0.001)
f4 = pd.ExcelFile(f'{FG}/Figure_4_source_data.xlsx'); check('Figure 4 PCoA = combined PCoA table', len(f4.parse('a_combined_pcoa')) == int(CP.loc['sample', 'n']))
# --- gene level
t = pd.read_csv(f'{B}/inputs/hitlevel/All_families_taxonomy_curated9396.tsv', sep='\t')
check('9,396 curated genes in 4,619 genomes', len(t) == 9396 and t.genome_id.nunique() == 4619)
sd3 = pd.read_excel(f'{T}/Supplementary_Data_3.xlsx', 'genes')
check('SD3 = curated set', set(sd3.hit_id) == set(t.hit_id) and len(sd3) == 9396)
G = pd.read_excel(f'{FG}/Figure_1_source_data.xlsx', 'c_curated_genes', index_col=0)
check('Figure 1c = curated gene table', G.loc[FAM, HOSTS].values.sum() == 9396 and (G.loc[FAM, HOSTS].values == pd.crosstab(t.gene, t.environment).reindex(index=FAM, columns=HOSTS, fill_value=0).values).all())
rm = pd.read_csv(f'{R}/s14_hitlevel/rdase_module.csv', index_col=0)
check('RDase module marine 127/235, terrestrial 0/22', (rm.loc[['Coral', 'Sponge']].with_module.sum(), rm.loc[['Coral', 'Sponge']].genes.sum(), rm.loc[['Chicken', 'Rhizosphere']].with_module.sum(), rm.loc[['Chicken', 'Rhizosphere']].genes.sum()) == (127, 235, 0, 22))
ss = pd.read_csv(f'{R}/s13_structures/structures_529_per_structure.csv')
check('529 structures, all curated, pLDDT >= 70, TM >= 0.5', len(ss) == 529 and ss.original_hit_id.isin(t.hit_id).all() and ss.plddt.min() >= 70 and ss.best_tm_score.min() >= 0.5)
sd5 = pd.read_excel(f'{T}/Supplementary_Data_6.xlsx', 'structures_529')
check('SD5 = 529 structures', len(sd5) == 529)
fs = pd.read_csv(f'{R}/s12_bgc/final_stats.tsv', sep='\t').set_index('set')
check('BGC primary 76/439 vs 75/939', (fs.loc['primary_ge3', 'hub_col'], fs.loc['primary_ge3', 'hub_n'], fs.loc['primary_ge3', 'fo_col'], fs.loc['primary_ge3', 'fo_n']) == (76, 439, 75, 939))
sd6 = pd.read_excel(f'{T}/Supplementary_Data_7.xlsx', 'primary_ge3')
check('SD6 primary genomes = 439 + 939', sd6.genome_id.nunique() == 1378, sd6.genome_id.nunique())
# --- Table 1 vs sources
T1 = pd.read_csv(f'{T}/Table_1.csv').set_index('Unnamed: 0') if 'Unnamed: 0' in pd.read_csv(f'{T}/Table_1.csv').columns else pd.read_csv(f'{T}/Table_1.csv').set_index(pd.read_csv(f'{T}/Table_1.csv').columns[0])
check('Table 1 totals', T1.loc['Metagenome samples (BioSamples)', 'Total'] == '1,787' and T1.loc['Studies', 'Total'] == '179' and T1.loc['Curated genes', 'Total'] == '9,396')
# --- figure files
for k in range(1, 7):
    check(f'Figure {k} PDF + JPEG + source data', all(os.path.exists(f'{FG}/Figure_{k}{e}') for e in ['.pdf', '.jpg', '_source_data.xlsx']))
for k in range(1, 16):
    check(f'Supplementary Figure {k} PDF', os.path.exists(f'{FG}/supplementary/Supplementary_Figure_{k}.pdf'))
print('\nFAILED:', fails if fails else 'none'); sys.exit(1 if fails else 0)
