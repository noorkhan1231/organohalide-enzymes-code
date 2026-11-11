"""Add round-2 panels to the formatted Supplementary Tables workbook (6 Nov 2026).
usage: python3 supp_tables_round2.py <Supplementary_Tables.xlsx in> <workbook out> <data archive folder>
Table S4 B-C: positions of unlabeled curated genes in the gene trees; Table S5 B-C: genes shared between host systems;
Table S10 E: chicken gut, leave-one-surroundings-project-out. Captions (row 2) updated."""
import sys, openpyxl, pandas as pd, numpy as np, math
from openpyxl.styles import Font, Alignment, PatternFill
inp, out, DEP = sys.argv[1:4]; DEP = DEP.rstrip('/'); RES = DEP + '/05_analysis_results'
wb = openpyxl.load_workbook(inp); NAVY = '1F3A5F'; hfill = PatternFill('solid', fgColor='EAF0F6')
def panel(ws, heading, desc, df):
    r = ws.max_row + 3; wid = 12
    ws.cell(r, 1, heading).font = Font(name='Calibri', size=12, bold=True, color=NAVY)
    c = ws.cell(r + 1, 1, desc); c.font = Font(name='Calibri', size=10, italic=True, color='444444'); c.alignment = Alignment(wrap_text=True, vertical='top')
    ws.merge_cells(start_row=r + 1, start_column=1, end_row=r + 1, end_column=wid); ws.row_dimensions[r + 1].height = 15 * max(1, math.ceil(len(desc) / (wid * 11)))
    for j, col in enumerate(df.columns, 1):
        x = ws.cell(r + 2, j, col); x.font = Font(name='Calibri', size=11, bold=True); x.fill = hfill
    for i, row in enumerate(df.itertuples(index=False), r + 3):
        for j, v in enumerate(row, 1):
            if isinstance(v, (float, np.floating)):
                if v != v: v = None
                elif np.isinf(v): v = 'inf'
                else: v = float(v)
            elif isinstance(v, (np.integer,)): v = int(v)
            elif isinstance(v, (np.bool_,)): v = bool(v)
            ws.cell(i, j, v)

# ---- S4
G = pd.read_csv(RES + '/s19_reference_free_clades/gene_positions.csv'); C = pd.read_csv(RES + '/s19_reference_free_clades/reference_free_clades.csv')
s = G.groupby(['family', 'unlabeled']).agg(curated_genes=('hit_id', 'size'), in_reference_free_clades=('in_ref_free_clade', 'sum'),
     median_patristic_distance_to_nearest_reference=('dist_nearest_ref', 'median'), median_identity_to_nearest_reference_pct=('identity', 'median')).reset_index()
s['label'] = np.where(s.unlabeled, 'unlabeled', 'labeled'); s['pct_in_reference_free_clades'] = (100 * s.in_reference_free_clades / s.curated_genes).round(1)
tot = G.groupby('unlabeled').agg(curated_genes=('hit_id', 'size'), in_reference_free_clades=('in_ref_free_clade', 'sum'),
     median_patristic_distance_to_nearest_reference=('dist_nearest_ref', 'median'), median_identity_to_nearest_reference_pct=('identity', 'median')).reset_index()
tot['family'] = 'all'; tot['label'] = np.where(tot.unlabeled, 'unlabeled', 'labeled'); tot['pct_in_reference_free_clades'] = (100 * tot.in_reference_free_clades / tot.curated_genes).round(1)
s = pd.concat([s, tot])[['family', 'label', 'curated_genes', 'in_reference_free_clades', 'pct_in_reference_free_clades', 'median_patristic_distance_to_nearest_reference', 'median_identity_to_nearest_reference_pct']].round(3)
cc = C.rename(columns={'curated_genes': 'curated_genes', 'labeled': 'labeled_genes', 'systems': 'host_systems'}).copy()
cc['unlabeled_genes'] = cc.curated_genes - cc.labeled_genes
cc = cc[['family', 'clade', 'curated_genes', 'labeled_genes', 'unlabeled_genes', 'host_systems'] + [c for c in cc.columns if c.startswith('n_')]].sort_values(['family', 'curated_genes'], ascending=[True, False])
ws = wb['Table S4']
panel(ws, 'B. Positions of labeled and unlabeled curated genes in the gene trees',
      'Unlabeled = curated genes without a subfamily label (unassigned, or families without subfamily resolution). Patristic distance (substitutions per site) from each curated gene to the nearest reference or outgroup tip of the same tree; reference-free clades = maximal clades holding at least ten curated genes and no reference or outgroup tip. Identity: DIAMOND identity to the closest characterized reference (Data S3).', s)
panel(ws, 'C. Reference-free clades with at least ten curated genes',
      'One row per clade: family, clade identifier, curated genes, labeled and unlabeled genes, host systems of the genes and counts per system. 44 clades hold no labeled gene (4,351 genes).', cc)
ws['A2'] = 'Panel A: composition of the ten gene trees. Panels B and C: positions of labeled and unlabeled curated genes in the trees and the clades that contain no characterized reference.'
# ---- S5
D = pd.read_excel(DEP + '/02_genes/Supplementary_Data_3.xlsx', sheet_name='genes')
from Bio import SeqIO
seq = {r.id: str(r.seq).rstrip('*') for r in SeqIO.parse(DEP + '/02_genes/curated_genes_protein.faa', 'fasta')}
D['seq'] = D.hit_id.map(seq)
ex = D.groupby('seq').agg(genes=('hit_id', 'size'), systems=('host_system', lambda v: '|'.join(sorted(set(v)))), family=('enzyme_family', 'first'), gene_ids=('hit_id', lambda v: ';'.join(v))).reset_index(drop=True)
ex = ex[ex.systems.str.contains(r'\|')]
P = pd.read_csv(RES + '/s18_shared_genes/cross_system_pairs_97.csv'); P['system_pair'] = ['|'.join(sorted([a, b])) for a, b in zip(P.system_a, P.system_b)]
summ = pd.concat([ex.groupby(['family', 'systems']).agg(sequences=('genes', 'size'), genes=('genes', 'sum')).reset_index().assign(level='identical protein').rename(columns={'systems': 'host_systems'}),
                  P.groupby(['family', 'system_pair']).agg(pairs=('gene_a', 'size'), genes=('gene_a', lambda v: len(set(v) | set(P.loc[v.index, 'gene_b'])))).reset_index().assign(level='>= 97% identity').rename(columns={'system_pair': 'host_systems'})])
summ = summ[['level', 'family', 'host_systems', 'sequences', 'pairs', 'genes']]
ws = wb['Table S5']
panel(ws, 'B. Curated genes shared between host systems',
      'No genome occurs in more than one host system. Identical protein: sequences carried by genes from two or more systems (22 genes). >= 97% identity: pairs of genes from different systems with at least 97% identity over the longer sequence in a global alignment (273 genes in 243 pairs; a lower bound, see STAR Methods).', summ)
panel(ws, 'C. Gene pairs from different host systems with at least 97% identity', 'One row per pair: family, gene identifiers, host systems and percent identity.', P[['family', 'gene_a', 'system_a', 'gene_b', 'system_b', 'identity']])
ws['A2'] = 'Panel A: counts of curated genes and of the genomes carrying them, with copy number and prevalence, for each family in each host system. Panels B and C: genes shared between host systems.'
# ---- S10
L = pd.read_csv(RES + '/s11_ambient/chicken_leave_one_project_out.csv')
ws = wb['Table S10']
panel(ws, 'E. Chicken gut versus surroundings, leaving out each surroundings project in turn',
      'project_left_out: surroundings BioProject removed (none = all ten projects). Mann–Whitney tests on study means as in panel C; P values are not adjusted. call_all_projects: call with all projects.', L.round(6))
ws['A2'] = 'Five panels. A: coverage of host and surrounding runs. B: null-gene recruitment. C: family abundance in host versus surrounding metagenomes with calls. D: all tests by scope and measure. E: chicken gut, leaving out each surroundings project in turn.'
wb.save(out); print('saved', out)
