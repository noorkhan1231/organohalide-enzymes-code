"""Step 15. Table 1, Supplementary Tables (one workbook) and Supplementary Data files from the final results.
usage: python3 s15_tables.py <FINAL dir> <old package dir (Bathy-RDases_ready_2026-10-21)> <out dir>
Only samples, runs and studies that are part of the final dataset appear; selection criteria are described in Methods.
"""
import sys, os, json, shutil, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stats_utils import HOSTS, FAM
B, OLD, O = sys.argv[1:4]; os.makedirs(O, exist_ok=True); R = f'{B}/results'; I = f'{B}/inputs'
S = pd.read_csv(f'{R}/s03_qc/samples.tsv', sep='\t'); S = S[S.set_main].copy()
rq = pd.read_csv(f'{R}/s03_qc/run_qc.tsv', sep='\t', low_memory=False)
rq = rq[(rq.exclude_code.isna() | rq.exclude_code.isin(['', 'LOWGE'])) & rq.sample_id.isin(S.sample_id)].copy()
CO = pd.read_csv(f'{R}/s08_map/sample_coordinates.tsv', sep='\t')
hl = json.load(open(f'{R}/s14_hitlevel/hitlevel_summary.json')); mt = json.load(open(f'{R}/s07_metatx/metatx_summary.json'))
F = pd.read_csv(f'{R}/s04_build/family_rpkg_main.csv')
W = pd.read_csv(f'{R}/s07_metatx/metatx_rpkm_per_run.csv')
AR = pd.read_csv(f'{R}/s11_ambient/ambient_runs_rpkg.csv')
assert len(S) == 1787 and S.sample_id.is_unique and set(rq.sample_id) == set(S.sample_id)

def med_iqr(v, scale=1, fmt='{:,.1f}'):
    v = np.asarray(v, float) / scale
    return f'{fmt.format(np.median(v))} ({fmt.format(np.percentile(v, 25))}–{fmt.format(np.percentile(v, 75))})'

# ================= Table 1 =================
amb = AR[AR.group != 'host']
AMBTYPE = {'Chicken': 'poultry litter and poultry-farm soil', 'Coral': 'seawater', 'Rhizosphere': 'bulk soil', 'Sponge': 'seawater'}
rows = []
def row(label, f, total=None): rows.append([label] + [f(h) for h in HOSTS] + [total])
sub = lambda h: S[S.host_system == h]
row('Metagenome samples (BioSamples)', lambda h: f'{len(sub(h)):,}', f'{len(S):,}')
row('Sequencing runs', lambda h: f'{(rq.host_system == h).sum():,}', f'{len(rq):,}')
row('Studies', lambda h: f'{sub(h).study.nunique()}', f'{S.study.nunique()}')
row('BioProjects', lambda h: f'{sub(h).BioProject.str.split(r"[;,]").explode().str.strip().nunique()}', f'{S.BioProject.str.split(r"[;,]").explode().str.strip().nunique()}')
cc = CO.set_index('sample_id').country
row('Countries or sea regions', lambda h: f'{cc[sub(h).sample_id].dropna().nunique()}', f'{cc.dropna().nunique()}')
row('Reads processed per sample, millions, median (IQR)', lambda h: med_iqr(sub(h).reads_processed, 1e6), med_iqr(S.reads_processed, 1e6))
row('Genome equivalents per sample, median (IQR)', lambda h: med_iqr(sub(h).GE, 1, '{:,.0f}'), med_iqr(S.GE, 1, '{:,.0f}'))
fr = F.set_index('sample_id')[FAM].sum(axis=1) > 0
row('Samples with reads on target genes, n (%)', lambda h: f'{fr[sub(h).sample_id].sum():,} ({100 * fr[sub(h).sample_id].mean():.0f})', f'{fr.sum():,} ({100 * fr.mean():.0f})')
row('Metatranscriptome runs (studies)', lambda h: f"{mt['runs_by_host'][h]} ({mt['studies_by_host'][h]})", f"{mt['n_runs']} ({W.study.nunique()})")
row('Surrounding metagenomes (projects)', lambda h: f"{amb[amb.system == h].run.nunique()} ({amb[amb.system == h].unit.nunique()})", f'{amb.run.nunique()} ({amb.unit.nunique()})')
row('Surrounding material', lambda h: AMBTYPE[h], '')
row('Genomes searched', lambda h: f"{hl['searched_by_host'][h]:,}", f"{hl['genomes_searched']:,}")
row('Genomes with curated genes', lambda h: f"{hl['carriers_by_host'][h]:,}", f"{hl['carrier_genomes']:,}")
G = pd.read_csv(f'{R}/s14_hitlevel/genes_per_family_host.csv', index_col=0)
row('Curated genes', lambda h: f'{G[h].sum():,}', f"{hl['curated_genes']:,}")
T1 = pd.DataFrame(rows, columns=['', *HOSTS, 'Total'])
with pd.ExcelWriter(f'{O}/Table_1.xlsx') as xw: T1.to_excel(xw, sheet_name='Table 1', index=False)
T1.to_csv(f'{O}/Table_1.csv', index=False)
print(T1.to_string(index=False))

# ================= Supplementary Tables =================
ST = {}; title = {}
o = pd.ExcelFile(f'{OLD}/Supplementary_Tables.xlsx')
for k, s in [('S1', 'S1_references_HMMs'), ('S2', 'S2_search_cutoffs'), ('S3', 'S3_gene_trees')]:
    ST[k] = o.parse(s)
_rf = pd.read_csv(f'{ROOT}/inputs/halocycdb_core_references.csv') if 'ROOT' in globals() else pd.read_csv('/home/claude/FINAL/inputs/halocycdb_core_references.csv')
_rf = _rf[(_rf.gene_family != 'HahyDase') & ~((_rf.gene_family == 'GSTDase') & (_rf.accession_number == 'MG878877.1'))]
assert len(_rf) == 217, len(_rf)
_rf = _rf[['gene_family', 'gene_name', 'accession_number', 'species', 'host', 'substrate', 'product', 'pfam', 'motif_conserved_seq', 'reference']].rename(columns={'gene_family': 'family', 'host': 'source_host'})
ST['S1'] = (ST['S1'], _rf)
title['S1'] = 'Characterised reference enzymes (217 accessions listed per family) and profile HMMs per enzyme family'
title['S2'] = 'Search cut-offs and genes per enzyme family and host system'
title['S3'] = 'Gene trees: tips, references, outgroups and curated genes'
# S4 studies
st = S.groupby('study').agg(host_system=('host_system', 'first'), BioProjects=('BioProject', lambda v: ', '.join(sorted(set(x.strip() for s_ in v for x in str(s_).replace(';', ',').split(','))))),
                            samples=('sample_id', 'size'), runs=('n_runs', 'sum'), median_GE=('GE', 'median'), median_reads_M=('reads_processed', lambda v: np.median(v) / 1e6)).reset_index()
cty = CO.groupby('study').country.agg(lambda v: '; '.join(sorted(set(v.dropna()))))
st['countries_or_sea_regions'] = st.study.map(cty)
ttl = rq.groupby('study').study_title.first(); st['study_title'] = st.study.map(ttl)
st['samples_with_target_reads'] = st.study.map(F.assign(r=fr.values).groupby('study').r.sum())
st = st.sort_values(['host_system', 'study']); ST['S4'] = st; title['S4'] = 'Studies in the final dataset'
# S5 family tests
K = pd.read_csv(f'{R}/s05_family/main_kruskal.csv'); LOO = pd.read_csv(f'{R}/s05_family/main_loo_kw.csv')
K['LOO_worst_P'] = [LOO[f].max() for f in K.family]; K['LOO_worst_study'] = [LOO.loc[LOO[f].idxmax(), 'dropped'] for f in K.family]
MT = pd.read_csv(f'{R}/s05_family/main_marine_vs_terrestrial.csv').add_prefix('marine_vs_terrestrial_').rename(columns={'marine_vs_terrestrial_family': 'family'})
ST['S5'] = K.merge(MT, on='family'); title['S5'] = 'Family abundance (RPKG) among host systems: Kruskal-Wallis at sample and study level, Dunn letters, leave-one-study-out, marine vs terrestrial'
# S6 community tests
V = pd.read_csv(f'{R}/s05_family/main_vegan_tests.csv'); V = V[V.test != 'host_restricted_within_block']
LG = pd.read_csv(f'{R}/s05_family/main_vegan_logGE_first.csv'); LG['test'] = 'log10(GE) entered first: ' + LG.term
SE = pd.read_csv(f'{R}/s05_family_sens/sensitivity_summary.csv')
IND = pd.read_csv(f'{R}/s05_family/main_indicator_families.csv'); DIS = pd.read_csv(f'{R}/s05_family/main_vegan_dispersion.csv')
ST['S6'] = (V, LG[['test', 'Df', 'R2', 'F', 'P']], SE, DIS, IND); title['S6'] = 'Composition of the ten families: PERMANOVA, PERMDISP, genome-equivalent covariate, genome-equivalent thresholds, indicator families'
RA = pd.read_csv(f'{R}/s05_family/main_relative_abundance.csv'); ST['S7'] = RA; title['S7'] = 'Relative abundance of the ten families within host systems (samples with target reads)'
PC = pd.read_csv(f'{R}/s14_hitlevel/prevalence_copy_number.csv'); ST['S8'] = PC; title['S8'] = 'Curated genes, carrier genomes, copy number and prevalence per family and host system'
SP = pd.read_csv(f'{R}/s06_subfamily/main_subfamily_permanova.csv'); SK = pd.read_csv(f'{R}/s06_subfamily/main_subfamily_kruskal.csv')
SI = pd.read_csv(f'{R}/s06_subfamily/main_subfamily_indicators_significant.csv')
HS = pd.read_csv(f'{R}/s16_indicator_support/indicator_host_vs_surroundings.csv'); HE = pd.read_csv(f'{R}/s16_indicator_support/indicator_expression.csv')
HL = pd.read_csv(f'{R}/s17_lineage_ranks/indicator_lineage_broader_ranks.csv')
ST['S9'] = (pd.read_csv(f'{R}/s06_subfamily/combined_permanova.csv'), SP, SI, SK, HS, HE, HL); title['S9'] = 'Subfamilies: composition by host system (PERMANOVA, all labelled subfamilies combined and per family), indicator subfamilies with genome-level lineage control, abundance tests, indicator subfamilies versus surrounding metagenomes and in metatranscriptomes, lineage control at GTDB family, class and phylum'
cov = pd.read_csv(f'{R}/s11_ambient/ambient_coverage.csv'); calls = pd.read_csv(f'{R}/s11_ambient/host_vs_ambient_calls.csv')
alltests = pd.read_csv(f'{R}/s11_ambient/host_vs_ambient_family.csv'); nul = pd.read_csv(f'{R}/s11_ambient/null_gene_recruitment.csv')
ST['S10'] = (cov, nul, calls, alltests); title['S10'] = 'Host-associated vs surrounding metagenomes: coverage, null-gene recruitment, calls and all tests'
MK = pd.read_csv(f'{R}/s07_metatx/metatx_family_tests.csv'); ST['S11'] = MK; title['S11'] = 'Metatranscriptomes: expression and detection of the ten families among host systems'
DV = pd.read_excel(f'{I}/hitlevel/Domain_vs_host_all_hits.xlsx', 'domain_vs_host'); DM = pd.read_excel(f'{I}/hitlevel/Domains_marine_terrestrial_and_pockets.xlsx', 'marine_vs_terrestrial')
DS = pd.read_excel(f'{I}/hitlevel/Subfamily_level_domain_and_genus_tests.xlsx', 'domain_by_subfamily')
DT = pd.read_csv(f'{R}/s14_hitlevel/deeptmhmm_summary.csv'); DTT = pd.read_csv(f'{I}/hitlevel/deeptmhmm_host_tests.tsv', sep='\t')
MB = pd.read_csv(f'{R}/s14_hitlevel/mobility_plasmid_contig.csv'); RM = pd.read_csv(f'{R}/s14_hitlevel/rdase_module.csv'); VM = pd.read_csv(f'{R}/s14_hitlevel/vanhase_module.csv')
ST['S12'] = (DV, DM, DS, RM, VM, DT, DTT, MB); title['S12'] = 'Domain architecture (Pfam-A), membrane topology (DeepTMHMM) and genomic context (geNomad) of curated genes'
ST['S13'] = (pd.read_csv(f'{R}/s12_bgc/SuppTable8_bgc_colocation_curated9396.csv'),); title['S13'] = 'Co-location of target genes with biosynthetic gene clusters (antiSMASH 6.1.1): comparisons and per-genus counts'
s13s = json.load(open(f'{R}/s13_structures/structure_summary_529.json'))
SS = pd.read_csv(f'{R}/s13_structures/structures_529_per_structure.csv')
sfam = SS.groupby('gene').agg(structures=('safe_id', 'size'), pLDDT_median=('plddt', 'median'), pLDDT_min=('plddt', 'min'), TM_median=('best_tm_score', 'median'),
                               TM_min=('best_tm_score', 'min'), identity_median=('percent_identity', 'median'), high_divergence_conserved_fold=('high_divergence_conserved_fold', 'sum')).reset_index()
ST['S14'] = (sfam, pd.read_csv(f'{R}/s13_structures/structure_permanova_529.csv'),
             pd.read_csv(f'{R}/s13_structures/pocket_length_adjusted_529.csv'),
             pd.read_csv(f'{R}/s13_structures/pocket_tool_agreement_529.csv'))
title['S14'] = 'Predicted structures (529): confidence, fold similarity, structural PERMANOVA by host system, pocket tests (raw, per 100 residues and adjusted for chain length) and agreement between the two pocket detectors'
ST['S15'] = (pd.read_csv(f'{R}/s10_geography/geo_mantel.csv'), pd.read_csv(f'{R}/s10_geography/geo_latitude.csv')); title['S15'] = 'Geography: Mantel tests (composition vs distance) and absolute latitude vs family abundance'
SW = pd.DataFrame([
    ('SRA Toolkit', '2.10.8', 'prefetch, fasterq-dump / fastq-dump'), ('fastp', '0.23.4', 'read trimming and filtering'), ('MicrobeCensus', '1.1.1', 'genome equivalents'),
    ('Salmon', '2.7.0', 'read mapping to target and null-gene indices'), ('HMMER', '3.3.2', 'hmmbuild, hmmsearch'), ('Clustal Omega', '1.2.4', 'alignments'), ('FastTreeMP', '2.2.0', 'gene trees (all ten families)'),
    ('CD-HIT', '4.8.1', 'sequence clustering for the structure sets'), ('Prodigal', '2.6.3', 'gene prediction'), ('gffread', '0.12.9', 'CDS extraction from genome annotations'), ('GTDB-Tk (GTDB R226)', '2.5.2', 'taxonomy'),
    ('DIAMOND', '2.1.8 (anvio-8 environment)', 'identity to references'), ('Pfam-A / hmmscan', 'Pfam 38.2 / HMMER', 'domain architecture'), ('DeepTMHMM', 'BioLib web server (run 23 Oct 2026)', 'topology'),
    ('geNomad', '1.11.2', 'plasmid/virus classification'), ('antiSMASH', '6.1.1', 'BGCs'), ('ColabFold', '1.6.3', 'structure prediction'),
    ('Foldseek', '10-941cd33', 'fold search'), ('fpocket', '4.2.3', 'pockets'), ('P2Rank', '2.5', 'pockets'),
    ('R', '4.3.3', 'statistics'), ('vegan', '2.6.4', 'PERMANOVA, PERMDISP, NMDS, Mantel'), ('Python', '3.13 (analysis), 3.10 (server)', 'analysis and table scripts'),
    ('pandas / NumPy / SciPy / statsmodels', '3.0.5 / 2.5.3 / 1.18.1 / 0.15.0', 'tables and tests'), ('matplotlib / GeoPandas', '3.11.2 / 0.14.4', 'figures, map (Natural Earth 1:110m)')],
    columns=['software', 'version', 'use'])
ST['S16'] = SW; title['S16'] = 'Software and versions'


# internal keys S1..S16 are the analysis order; the published numbering follows first citation in the manuscript
TMAP = {4: 1, 1: 2, 2: 3, 3: 4, 8: 5, 6: 6, 7: 7, 5: 8, 15: 9, 10: 10, 11: 11, 9: 12, 12: 13, 14: 14, 13: 15, 16: 16}
NEW = sorted(ST, key=lambda k: TMAP[int(k[1:])]); nn = {k: f'S{TMAP[int(k[1:])]}' for k in ST}
with pd.ExcelWriter(f'{O}/Supplementary_Tables.xlsx') as xw:
    pd.DataFrame({'Sheet': [f'Table {nn[k]}' for k in NEW], 'Title': [title[k] for k in NEW]}).to_excel(xw, sheet_name='Contents', index=False)
    for k in NEW:
        v = ST[k]; parts = v if isinstance(v, tuple) else (v,); r0 = 2
        pd.DataFrame({f'Supplementary Table {nn[k][1:]}. {title[k]}': []}).to_excel(xw, sheet_name=f'Table {nn[k]}', index=False)
        for p in parts:
            p.to_excel(xw, sheet_name=f'Table {nn[k]}', index=False, startrow=r0); r0 += len(p) + 3
print('Supplementary Tables:', list(ST))

# ================= Supplementary Data =================
# SD1 samples and runs
sd1s = S[['sample_id', 'host_system', 'study', 'BioProject', 'runs', 'n_runs', 'reads_processed', 'Gb', 'GE', 'platform']].merge(
    CO[['sample_id', 'lat', 'lon', 'country']], on='sample_id', how='left')
sd1s['reads_on_target_genes'] = sd1s.sample_id.map(pd.read_csv(f'{R}/s04_build/family_reads_main.csv').set_index('sample_id')[FAM].sum(axis=1))
cols_r = ['run', 'sample_id', 'host_system', 'study', 'BioProject', 'sample_accession', 'scientific_name', 'host', 'isolation_source', 'library_strategy', 'library_source',
          'instrument_platform', 'reads_processed', 'GE', 'fastp_reads_in', 'fastp_reads_out', 'q30_after', 'collection_date', 'study_title']
sd1r = rq[[c for c in cols_r if c in rq.columns]]
with pd.ExcelWriter(f'{O}/Supplementary_Data_1.xlsx') as xw:
    pd.DataFrame({'Supplementary Data 1. Metagenome samples and runs': [f'samples: {len(sd1s):,} BioSamples (runs of one BioSample pooled); runs: {len(sd1r):,} sequencing runs. lat and lon give the sample location identified from the sample records (11 samples have none). country is blank where the record names no country or sea region (85 samples, mostly coral). library_strategy and library_source are the labels in the ENA record; library_source is blank where the record gives none. Fifty-five sponge runs carry the labels AMPLICON (PRJNA1280067, 32 runs) or OTHER (PRJNA1256488, 22 runs; PRJNA1072552, 1 run) but are described by their NCBI project records as shotgun metagenomes of sponge tissue and were retained.']}).to_excel(xw, sheet_name='README', index=False)
    sd1s.to_excel(xw, sheet_name='samples', index=False); sd1r.to_excel(xw, sheet_name='runs', index=False)
# SD2 gene trees (unchanged)
shutil.copy(f'{OLD}/Supplementary_Data/Supplementary_Data_2.zip', f'{O}/Supplementary_Data_2.zip')
# SD3 curated genes
o3 = pd.read_excel(f'{OLD}/Supplementary_Data/Supplementary_Data_3.xlsx', 'hits')
A = pd.read_excel(f'{I}/hitlevel/Domain_architecture_per_hit.xlsx')[['hit_id', 'subfamily_call', 'arch']].rename(columns={'subfamily_call': 'subfamily', 'arch': 'Pfam_architecture'})
D = pd.read_csv(f'{I}/hitlevel/deeptmhmm_per_hit_annotated.tsv', sep='\t')[['hit_id', 'dtm_type', 'len']].rename(columns={'dtm_type': 'DeepTMHMM_topology', 'len': 'length_aa'})
sd3 = o3.merge(A, on='hit_id', how='left').merge(D, on='hit_id', how='left')
assert len(sd3) == 9396 and sd3.subfamily.notna().all() and sd3.DeepTMHMM_topology.notna().all()
with pd.ExcelWriter(f'{O}/Supplementary_Data_3.xlsx') as xw:
    pd.DataFrame({'Supplementary Data 3. Curated genes of the ten enzyme families': ['One row per curated gene (9,396 genes in 4,619 genomes): source genome, host system, GTDB lineage, identity to the closest characterised reference, subfamily label (unassigned = no label; no_motif_signal = family without subfamily resolution), Pfam-A architecture (hmmscan --cut_ga, i-Evalue <= 1e-5) and DeepTMHMM topology (SP signal peptide, TM transmembrane, GLOB globular).']}).to_excel(xw, sheet_name='README', index=False)
    sd3.to_excel(xw, sheet_name='genes', index=False)
# SD4 abundance
with pd.ExcelWriter(f'{O}/Supplementary_Data_4.xlsx') as xw:
    pd.DataFrame({'Supplementary Data 4. Abundance of the ten families and their subfamilies per sample': ['RPKG = reads per kilobase of target per genome equivalent (all 9,187 mapping targets; runs of one BioSample pooled). family_reads: reads assigned to each family. subfamily_RPKG: long format, samples with reads only.']}).to_excel(xw, sheet_name='README', index=False)
    F.to_excel(xw, sheet_name='family_RPKG', index=False)
    pd.read_csv(f'{R}/s04_build/family_reads_main.csv').to_excel(xw, sheet_name='family_reads', index=False)
    pd.read_csv(f'{R}/s04_build/subfamily_rpkg_long_main.csv').to_excel(xw, sheet_name='subfamily_RPKG', index=False)
# SD5 structures (529)
x5 = pd.ExcelFile(f'{OLD}/Supplementary_Data/Supplementary_Data_5.xlsx'); keep = set(SS.safe_id)
with pd.ExcelWriter(f'{O}/Supplementary_Data_6.xlsx') as xw:
    pd.DataFrame({'Supplementary Data 6. Predicted structures of curated genes (529)': ['Per-structure confidence (pLDDT), pocket metrics (fpocket, P2Rank), Foldseek best matches and Pfam domains. Tests by host system are in Supplementary Table 14. In p2rank_pockets, P2Rank scores are blank for the one structure (VanHase_V_IPO_Coral_13) in which P2Rank predicted no pocket.']}).to_excel(xw, sheet_name='README', index=False)
    SS.to_excel(xw, sheet_name='structures_529', index=False)
    for sh in ['four_groups_plddt_pockets', 'GSTDase_plddt_pockets', 'MetDase_NHFeHase_plddt_pockets', 'fold_search_flagship', 'fold_search_indicator', 'fold_search_GSTDase',
               'fold_search_MetDase_NHFeHase', 'pfam_domains', 'p2rank_pockets']:
        d = x5.parse(sh); idc = 'safe_id' if 'safe_id' in d.columns else d.columns[0]
        d = d[d[idc].isin(keep)]; d.to_excel(xw, sheet_name=sh, index=False)
# SD6 BGC regions
with pd.ExcelWriter(f'{O}/Supplementary_Data_7.xlsx') as xw:
    pd.DataFrame({'Supplementary Data 7. Biosynthetic gene clusters in genomes analysed for co-location': ['antiSMASH 6.1.1 regions with product class, nearest target gene and relationship (within, flanking <= 10 kb, distant). primary: hub genera with labels in >= 3 families; sensitivity: >= 2.']}).to_excel(xw, sheet_name='README', index=False)
    for s_ in ['primary_ge3', 'sensitivity_ge2']:
        pd.read_csv(f'{I}/bgc_genomic_hubs/{s_}/bgc_per_genome.tsv', sep='\t').to_excel(xw, sheet_name=s_, index=False)
# SD7 metatranscriptomes and surroundings
with pd.ExcelWriter(f'{O}/Supplementary_Data_5.xlsx') as xw:
    pd.DataFrame({'Supplementary Data 5. Metatranscriptomes and surrounding metagenomes': ['metatranscriptomes: reads per kilobase per million reads processed per family (first 30 million reads). surroundings: RPKG per family and null-gene RPKG of each surrounding metagenome; host_null: null-gene RPKG of host runs.']}).to_excel(xw, sheet_name='README', index=False)
    W.to_excel(xw, sheet_name='metatranscriptomes', index=False)
    AR[AR.group != 'host'].drop(columns=['ds', 'sample_id']).rename(columns={'unit': 'BioProject'}).to_excel(xw, sheet_name='surroundings', index=False)
    AR[(AR.group == 'host') & AR.null_RPKG.notna()][['run', 'system', 'unit', 'GE', 'nproc_full', 'nproc_sub', 'null_rpk_sub', 'null_RPKG']].rename(columns={'unit': 'study'}).to_excel(xw, sheet_name='host_null', index=False)
print('Supplementary Data written to', O)
