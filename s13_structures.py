"""Step 12. Predicted-structure numbers on the curated gene set (E9: 537 -> 529 structures).
usage: python3 s13_structures.py <Supplementary_Data_5.xlsx> <All_families_taxonomy_curated9396.tsv> <struct_permanova_curated.csv> <out_dir>
Structures whose source hit is not among the 9,396 curated genes are removed (they passed the HMM thresholds but were
removed later by phylogenetic curation). Recomputes: counts per family/set, pLDDT summary, Foldseek best-match TM-score
summary, high-divergence/conserved-fold and candidate-novel-fold counts, and fpocket/P2Rank Kruskal-Wallis tests by host
within each subfamily group (BH within metric). Structural PERMANOVA (TM-score dissimilarity) on the curated set comes from
struct_permanova_curated.csv (vegan adonis2, 999 permutations; computed with build_dist.py + perm.R).
"""
import sys, os, numpy as np, pandas as pd
from scipy.stats import kruskal
from statsmodels.stats.multitest import multipletests
sd5, taxf, permf, O = sys.argv[1:5]; os.makedirs(O, exist_ok=True)
x = pd.ExcelFile(sd5)
cur = set(pd.read_csv(taxf, sep='\t', usecols=['hit_id']).hit_id)
fs = pd.concat([x.parse(s).assign(sheet=s) for s in ['fold_search_flagship', 'fold_search_indicator', 'fold_search_GSTDase', 'fold_search_MetDase_NHFeHase']], ignore_index=True)
fs['curated'] = fs.original_hit_id.isin(cur)
rem = fs[~fs.curated]
print('structures', len(fs), '| not in curated set:', len(rem)); print(rem[['safe_id', 'gene', 'subfamily', 'environment', 'original_hit_id', 'sheet']].to_string(index=False))
rem.to_csv(f'{O}/structures_removed_E9.csv', index=False)
keep = set(fs[fs.curated].safe_id)
pl = pd.concat([x.parse('four_groups_plddt_pockets').rename(columns={'mean_plddt': 'plddt', 'n_pockets_detected': 'n_pockets'})[['safe_id', 'gene', 'plddt', 'n_pockets', 'n_residues']],
                x.parse('GSTDase_plddt_pockets')[['safe_id', 'gene', 'plddt', 'n_pockets', 'n_residues']],
                x.parse('MetDase_NHFeHase_plddt_pockets')[['safe_id', 'gene', 'plddt', 'n_pockets', 'n_residues']]], ignore_index=True)
pl = pl[pl.safe_id.isin(keep)]
f2 = fs[fs.curated]
S = dict(structures=len(f2), by_family=f2.gene.value_counts().to_dict(),
         plddt_min=round(pl.plddt.min(), 1), plddt_median=round(pl.plddt.median(), 1), plddt_median_by_family=pl.groupby('gene').plddt.median().round(1).to_dict(),
         plddt_ge70=int((pl.plddt >= 70).sum()), plddt_ge90=int((pl.plddt >= 90).sum()),
         tm_min=round(f2.best_tm_score.min(), 4), tm_median=round(f2.best_tm_score.median(), 4), tm_median_by_family=f2.groupby('gene').best_tm_score.median().round(3).to_dict(),
         tm_below_0_5=int((f2.best_tm_score < 0.5).sum()), candidate_novel_fold=int(f2.candidate_novel_fold.sum()),
         high_divergence_conserved_fold=int(f2.high_divergence_conserved_fold.sum()),
         identity_median=round(f2.percent_identity.median(), 1), identity_below_30=int((f2.percent_identity < 30).sum()),
         best_database=f2.best_database.value_counts().to_dict())
pd.Series(S).to_json(f'{O}/structure_summary_529.json', indent=1)
per = f2[['safe_id', 'gene', 'subfamily', 'environment', 'original_hit_id', 'best_tm_score', 'percent_identity', 'best_database',
          'candidate_novel_fold', 'high_divergence_conserved_fold']].merge(pl[['safe_id', 'plddt', 'n_residues']], on='safe_id', how='left')
assert per.plddt.notna().all() and len(per) == len(f2), 'every curated structure needs a pLDDT value'
per.to_csv(f'{O}/structures_529_per_structure.csv', index=False)
for k, v in S.items(): print(k, v)
# pocket tests by host within subfamily group, curated structures only (fpocket counts per 100 residues and P2Rank metrics)
p2 = x.parse('p2rank_pockets'); p2 = p2[p2.safe_id.isin(keep)]
rows = []
for (fam, sub), g in p2.groupby(['Enzyme family', 'Subfamily']):
    if g['Host system'].nunique() < 2 or g.groupby('Host system').size().min() < 2: continue
    for m in ['fpocket pockets per 100 residues', 'P2Rank pockets per 100 residues', 'P2Rank mean score, top 3 pockets', 'Chain length (residues)']:
        vals = [v[m].dropna().values for _, v in g.groupby('Host system')]
        if sum(len(v) for v in vals) < 4 or np.ptp(np.concatenate(vals)) == 0: continue
        h, p = kruskal(*vals); rows.append(dict(family=fam, subfamily=sub, metric=m, n=len(g), hosts=g['Host system'].nunique(), H=h, P=p))
T = pd.DataFrame(rows)
T['FDR'] = np.nan
for m in T.metric.unique():
    k = T.metric == m; T.loc[k, 'FDR'] = multipletests(T.loc[k, 'P'], method='fdr_bh')[1]
T.to_csv(f'{O}/pocket_tests_529.csv', index=False)
print('\npocket tests with FDR < 0.05:'); print(T[T.FDR < 0.05].round(4).to_string(index=False))
P = pd.read_csv(permf); P['FDR'] = np.nan; k = P.P.notna(); P.loc[k, 'FDR'] = multipletests(P.loc[k, 'P'], method='fdr_bh')[1]
k2 = P.P_given_length.notna(); P.loc[k2, 'FDR_given_length'] = multipletests(P.loc[k2, 'P_given_length'], method='fdr_bh')[1]
P.to_csv(f'{O}/structure_permanova_529.csv', index=False)
print('\nstructural PERMANOVA (curated):'); print(P[P.P.notna()].round(4).to_string(index=False))
