"""Chicken gut vs surroundings: leave-one-surroundings-project-out (Prof. Zhang round 2; 6 Nov 2026).
Same test as s11_ambient.py (two-sided Mann-Whitney on study means; RPKG on all main-set host samples, null-normalized on runs with a null value)."""
import sys, numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
B = sys.argv[1]; R = f'{B}/results'
FAM = ['RDase', 'HyDase', 'GSTDase', 'OxDase', 'DehyDase', 'MetDase', 'FlaHase', 'VanHase', 'NHFeHase', 'SAMHase']
A = pd.read_csv(f'{R}/s11_ambient/ambient_runs_rpkg.csv'); A = A[A.system == 'Chicken']
F = pd.read_csv(f'{R}/s04_build/family_rpkg_main.csv'); H = F[F.host_system == 'Chicken']
C = pd.read_csv(f'{R}/s11_ambient/host_vs_ambient_calls.csv'); C = C[C.system == 'Chicken'].set_index('set')
amb = A[A.group != 'host']; projects = sorted(amb.unit.unique()); rows = []
for f in FAM:
    h = H.groupby('study')[f].mean(); a = amb.groupby('unit')[f].mean()
    x = A[A.null_RPKG > 0].copy(); x['n'] = x[f] / x.null_RPKG
    hn = x[x.group == 'host'].groupby('unit').n.mean(); an = x[x.group != 'host'].groupby('unit').n.mean()
    for u in ['none'] + projects:
        aa = a if u == 'none' else a.drop(u); ann = an if u == 'none' else an.drop(u)
        if (h.sum() + aa.sum()) == 0: continue
        rows.append(dict(family=f, call_all_projects=C.loc[f, 'call'], project_left_out=u, n_host_studies=len(h), n_surroundings_projects=len(aa),
                         ratio_RPKG=h.mean() / aa.mean() if aa.mean() > 0 else np.inf, P_RPKG=mannwhitneyu(h, aa, alternative='two-sided').pvalue,
                         ratio_null_normalized=hn.mean() / ann.mean() if ann.mean() > 0 else np.inf, P_null_normalized=mannwhitneyu(hn, ann, alternative='two-sided').pvalue))
T = pd.DataFrame(rows)
T['same_direction_both'] = [(r.ratio_RPKG < 1) == (r.ratio_null_normalized < 1) for r in T.itertuples()]
T.to_csv(f'{R}/s11_ambient/chicken_leave_one_project_out.csv', index=False)
S = T[T.project_left_out != 'none'].groupby(['family', 'call_all_projects']).agg(
    runs=('project_left_out', 'size'), max_P_RPKG=('P_RPKG', 'max'), max_P_null_normalized=('P_null_normalized', 'max'),
    host_lower_in_all_RPKG=('ratio_RPKG', lambda v: bool((v < 1).all())), host_lower_in_all_null=('ratio_null_normalized', lambda v: bool((v < 1).all()))).reset_index()
S.to_csv(f'{R}/s11_ambient/chicken_leave_one_project_out_summary.csv', index=False)
print(S.round(5).to_string(index=False))
