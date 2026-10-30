import numpy as np, pandas as pd, itertools, subprocess, os
from scipy.stats import kruskal, rankdata, norm, spearmanr
from statsmodels.stats.multitest import multipletests

HOSTS = ['Chicken', 'Coral', 'Rhizosphere', 'Sponge']
FAM = ['RDase', 'HyDase', 'GSTDase', 'OxDase', 'DehyDase', 'MetDase', 'FlaHase', 'VanHase', 'NHFeHase', 'SAMHase']

def bh(p):
    p = np.asarray(p, float); out = np.full(len(p), np.nan); ok = ~np.isnan(p)
    if ok.any(): out[ok] = multipletests(p[ok], method='fdr_bh')[1]
    return out

def kw(values, groups):
    vals = [values[groups == g] for g in HOSTS if (groups == g).any()]
    if sum(len(v) > 0 for v in vals) < 2 or np.ptp(np.concatenate(vals)) == 0: return np.nan, np.nan
    r = kruskal(*vals); return r.statistic, r.pvalue

def dunn(values, groups):
    """Dunn's test with tie correction, BH-adjusted; returns list of (g1, g2, z, p, p_adj)."""
    x = np.asarray(values, float); g = np.asarray(groups); N = len(x); r = rankdata(x)
    _, cnt = np.unique(x, return_counts=True); tie = (cnt ** 3 - cnt).sum() / (12 * (N - 1))
    out = []
    for a, b in itertools.combinations([h for h in HOSTS if (g == h).any()], 2):
        na, nb = (g == a).sum(), (g == b).sum()
        se = np.sqrt((N * (N + 1) / 12 - tie) * (1 / na + 1 / nb))
        z = (r[g == a].mean() - r[g == b].mean()) / se if se > 0 else 0.0
        out.append([a, b, z, 2 * norm.sf(abs(z))])
    if out:
        adj = bh([o[3] for o in out]); [o.append(q) for o, q in zip(out, adj)]
    return out

def cld(groups, pairs, means, alpha=0.05):
    """Compact letter display from pairwise adjusted P (insert-and-absorb), letters ordered by mean (a = highest)."""
    gs = sorted(groups, key=lambda h: -means[h]); sig = {(a, b): q < alpha for a, b, *_, q in pairs}
    sig.update({(b, a): v for (a, b), v in list(sig.items())})
    sets = [set(gs)]
    for a, b in itertools.combinations(gs, 2):
        if sig.get((a, b)):
            new = []
            for s in sets:
                if a in s and b in s: new += [s - {a}, s - {b}]
                else: new.append(s)
            new = [s for s in new if not any(s < t for t in new)]
            sets = list({frozenset(s): s for s in new}.values())
    sets = sorted(sets, key=lambda s: min(gs.index(h) for h in s))
    letters = {h: '' for h in gs}
    for i, s in enumerate(sets):
        for h in s: letters[h] += 'abcdefghij'[i]
    return letters

def shannon(M):
    P = M / M.sum(axis=1, keepdims=True)
    with np.errstate(divide='ignore', invalid='ignore'): return -np.nansum(np.where(P > 0, P * np.log(P), 0), axis=1)

def run_vegan(df, id_col, group_col, block_col, desc_col, prefix, k):
    f = prefix + '_input.csv'; df.to_csv(f, index=False)
    cmd = ['Rscript', os.path.join(os.path.dirname(__file__), 'vegan_tests.R'), f, id_col, group_col, block_col, desc_col, prefix, str(k)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode: raise RuntimeError(r.stderr[-2000:])
    return pd.read_csv(prefix + '_tests.csv')

def cooccur(M, cols, min_frac=0.0):
    rows = []
    for a, b in itertools.combinations(range(M.shape[1]), 2):
        x, y = M[:, a], M[:, b]
        if np.ptp(x) == 0 or np.ptp(y) == 0: continue
        rho, p = spearmanr(x, y); rows.append((cols[a], cols[b], rho, p))
    d = pd.DataFrame(rows, columns=['a', 'b', 'rho', 'p']); d['fdr'] = bh(d.p.values) if len(d) else []
    return d
