"""Step 8. Family-level statistics on sample-level RPKG (all targets).
usage: python3 s05_family.py <family_rpkg_<set>.csv> <outdir> [label]
Sample = pooled BioSample; study = BioProject or merged study group (study_groups.tsv).
Tests: Kruskal-Wallis per family (sample level and study level, BH FDR), Dunn post hoc with letters, leave-one-study-out,
Shannon diversity, relative composition, indicator families (multipatt), within-host co-occurrence, PERMANOVA + PERMDISP
(sample level, nested in study, permutations restricted within study, study level), PERMANOVA with log10 genome equivalents
entered first, and the marine (Coral + Sponge) vs terrestrial (Chicken + Rhizosphere) contrast."""
import sys, os, json, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stats_utils import *
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lib'))
from multipatt import multipatt
fin, out = sys.argv[1:3]; label = sys.argv[3] if len(sys.argv) > 3 else 'main'; desc = '-'
os.makedirs(out, exist_ok=True); P = f'{out}/{label}'
d = pd.read_csv(fin)
tot = d[FAM].sum(axis=1); zero = d[tot == 0]          # samples without any target read stay in abundance tests (RPKG = 0);
# they are left out only where a composition is needed (Shannon, relative abundance, Bray-Curtis tests)
R = {'n_quantified': int(len(d)), 'n_zero': int(len(zero)), 'zero_by_host': zero.host_system.value_counts().to_dict(),
     'n': int(len(d)), 'n_with_reads': int((tot > 0).sum()), 'n_by_host': d.host_system.value_counts().reindex(HOSTS).to_dict(),
     'bioprojects_by_host': d.groupby('host_system').BioProject.nunique().reindex(HOSTS).to_dict(),
     'studies_by_host': d.groupby('host_system').study.nunique().reindex(HOSTS).to_dict()}
g = d.host_system.values
st = d.groupby(['host_system', 'study'])[FAM].mean().reset_index()
rows, posthoc = [], []
for f in FAM:
    H, p = kw(d[f].values, g); Hs, ps = kw(st[f].values, st.host_system.values)
    rows.append(dict(family=f, H=H, P=p, H_study=Hs, P_study=ps, **{f'mean_{h}': d.loc[g == h, f].mean() for h in HOSTS},
                     **{f'median_{h}': d.loc[g == h, f].median() for h in HOSTS}, **{f'prev_{h}': (d.loc[g == h, f] > 0).mean() for h in HOSTS}))
    pr = dunn(d[f].values, g); rk = pd.Series(d[f]).rank().values
    let = cld(HOSTS, pr, {h: rk[g == h].mean() for h in HOSTS})
    for a, b, z, pp, q in pr: posthoc.append(dict(family=f, a=a, b=b, z=z, P=pp, FDR=q))
    rows[-1].update({f'letter_{h}': let[h] for h in HOSTS})
K = pd.DataFrame(rows); K['FDR'] = bh(K.P.values); K['FDR_study'] = bh(K.P_study.values)
K.to_csv(f'{P}_kruskal.csv', index=False); pd.DataFrame(posthoc).to_csv(f'{P}_dunn.csv', index=False)
R['kw_n_sig'] = int((K.FDR < 0.05).sum()); R['kw_max_FDR'] = float(K.FDR.max()); R['kw_study_n_sig'] = int((K.FDR_study < 0.05).sum())
R['kw_study_sig'] = K.loc[K.FDR_study < 0.05, 'family'].tolist()
# leave-one-study-out: worst-case KW P per family
loo = []
for b in sorted(d.study.unique()):
    k = d.study != b
    if d.loc[k].groupby('host_system').study.nunique().reindex(HOSTS).fillna(0).min() < 1: continue
    loo.append(dict(dropped=b, host=d.loc[~k, 'host_system'].iat[0], **{f: kw(d.loc[k, f].values, d.loc[k, 'host_system'].values)[1] for f in FAM}))
L = pd.DataFrame(loo); L.to_csv(f'{P}_loo_kw.csv', index=False)
R['loo_worst_P'] = {f: float(L[f].max()) for f in FAM}; R['loo_worst_project'] = {f: L.loc[L[f].idxmax(), 'dropped'] for f in FAM}
# diversity and composition within hosts
dn = d[tot > 0].reset_index(drop=True); gn = dn.host_system.values
dn['shannon'] = shannon(dn[FAM].values); Hsh, psh = kw(dn.shannon.values, gn)
R['shannon_mean'] = dn.groupby('host_system').shannon.mean().reindex(HOSTS).round(3).to_dict(); R['shannon_H'] = Hsh; R['shannon_P'] = psh
rel = dn[FAM].div(dn[FAM].sum(axis=1), axis=0) * 100; rel['host_system'] = gn
top = dn[FAM].idxmax(axis=1)
RA = []
for h in HOSTS:
    x = rel[rel.host_system == h][FAM]
    for f in FAM: RA.append(dict(host=h, family=f, n=len(x), mean=x[f].mean(), median=x[f].median(), q1=x[f].quantile(.25), q3=x[f].quantile(.75), max=x[f].max(), pct_top=(top[gn == h] == f).mean() * 100))
pd.DataFrame(RA).to_csv(f'{P}_relative_abundance.csv', index=False)
# indicator families (all community samples)
lev, comb, stat, pv = multipatt(d[FAM].values, list(g), 999, 42)
I = pd.DataFrame(dict(family=FAM, combination=['+'.join(c) for c in comb], stat=stat, P=pv)); I['FDR'] = bh(I.P.values)
I.to_csv(f'{P}_indicator_families.csv', index=False); R['indicator_families'] = I[I.FDR <= 0.05][['family', 'combination', 'stat']].round(3).values.tolist()
# co-occurrence within host
C = []
for h in HOSTS:
    x = d[g == h]; c = cooccur(x[FAM].values, FAM); c['host'] = h; C.append(c)
C = pd.concat(C); C.to_csv(f'{P}_cooccurrence.csv', index=False)
sigc = C[C.fdr < 0.05]; R['cooc_sig_by_host'] = sigc.host.value_counts().reindex(HOSTS).fillna(0).astype(int).to_dict()
R['cooc_negative'] = sigc[sigc.rho < 0][['host', 'a', 'b', 'rho']].round(3).values.tolist()
# vegan
vd = d[['sample_id', 'host_system', 'study'] + FAM]
V = run_vegan(vd, 'sample_id', 'host_system', 'study', desc, P + '_vegan', int(os.environ.get('NMDS_K', '2')))
R['vegan'] = V.round(5).to_dict('records')
disp = pd.read_csv(P + '_vegan_dispersion.csv'); R['dispersion'] = dict(zip(disp.group, disp.mean_dist.round(3)))
# PERMANOVA with log10 genome equivalents entered before host (sequential terms)
import subprocess
d[['sample_id', 'host_system', 'GE'] + FAM].to_csv(P + '_ge_input.csv', index=False)
rc = '''suppressMessages(library(vegan)); d <- read.csv(commandArgs(TRUE)[1]); X <- d[, 4:13]; k <- rowSums(X) > 0; d <- d[k, ]; X <- X[k, ]
set.seed(42); p <- adonis2(vegdist(X, "bray") ~ log10(GE) + host_system, data = d, permutations = 999, by = "terms")
write.csv(data.frame(term = rownames(p), Df = p$Df, R2 = p$R2, F = p$F, P = p$`Pr(>F)`), commandArgs(TRUE)[2], row.names = FALSE)'''
open(P + '_ge.R', 'w').write(rc)
r = subprocess.run(['Rscript', P + '_ge.R', P + '_ge_input.csv', P + '_vegan_logGE_first.csv'], capture_output=True, text=True)
if r.returncode: raise RuntimeError(r.stderr[-1500:])
R['permanova_logGE_first'] = pd.read_csv(P + '_vegan_logGE_first.csv').round(5).to_dict('records')
# marine (Coral + Sponge) vs terrestrial (Chicken + Rhizosphere): Mann-Whitney per family, sample and study level
from scipy.stats import mannwhitneyu
mar = d.host_system.isin(['Coral', 'Sponge']).values; smar = st.host_system.isin(['Coral', 'Sponge']).values
MT = []
for f in FAM:
    u, p = mannwhitneyu(d.loc[mar, f], d.loc[~mar, f]); us, ps = mannwhitneyu(st.loc[smar, f], st.loc[~smar, f])
    MT.append(dict(family=f, mean_marine=d.loc[mar, f].mean(), mean_terrestrial=d.loc[~mar, f].mean(), median_marine=d.loc[mar, f].median(), median_terrestrial=d.loc[~mar, f].median(),
                   prev_marine=(d.loc[mar, f] > 0).mean(), prev_terrestrial=(d.loc[~mar, f] > 0).mean(), U=u, P=p, U_study=us, P_study=ps))
MT = pd.DataFrame(MT); MT['FDR'] = bh(MT.P.values); MT['FDR_study'] = bh(MT.P_study.values); MT.to_csv(f'{P}_marine_vs_terrestrial.csv', index=False)
R['marine_vs_terrestrial_sig'] = MT.loc[MT.FDR < 0.05, 'family'].tolist(); R['marine_vs_terrestrial_sig_study'] = MT.loc[MT.FDR_study < 0.05, 'family'].tolist()
json.dump(R, open(f'{P}_summary.json', 'w'), indent=1, default=float)
print(json.dumps(R, indent=1, default=float)[:6000])
print(K[['family', 'H', 'FDR', 'H_study', 'FDR_study'] + [f'letter_{h}' for h in HOSTS]].round(5).to_string(index=False))
