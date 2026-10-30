# Python re-implementation of indicspecies::multipatt(func="IndVal.g", duleg=FALSE), all group combinations except the full set.
import numpy as np, itertools
def combos(G):
    out=[]
    for r in range(1,G+1):
        for c in itertools.combinations(range(G),r): out.append(c)
    return out   # order matches multipatt index (singles, pairs, triples ...)
def best_stat(X, gi, G, C):
    # X: sites x species (abundance), gi: group index per site
    M=np.vstack([X[gi==k].mean(0) for k in range(G)])          # G x S mean abundance
    P=np.vstack([(X[gi==k]>0).sum(0) for k in range(G)])       # G x S presences
    N=np.array([(gi==k).sum() for k in range(G)])[:,None]
    tot=M.sum(0); tot[tot==0]=np.nan
    st=[]
    for c in C:
        A=M[list(c)].sum(0)/tot; B=P[list(c)].sum(0)/N[list(c)].sum()
        st.append(np.sqrt(A*B))
    st=np.vstack(st)                                           # combos x S
    return np.nanmax(st,0), np.nanargmax(np.nan_to_num(st,nan=-1),0)
def multipatt(X, groups, nperm=999, seed=42):
    lev=sorted(set(groups)); gi=np.array([lev.index(g) for g in groups]); G=len(lev); C=combos(G)
    obs,idx=best_stat(X,gi,G,C)
    rng=np.random.default_rng(seed); cnt=np.zeros(X.shape[1])
    for _ in range(nperm):
        p,_=best_stat(X,rng.permutation(gi),G,C); cnt+= (p>=obs-1e-12)
    pval=(cnt+1)/(nperm+1)
    pval=np.where((idx==len(C)-1)|np.isnan(obs), np.nan, pval)   # best = all groups -> not tested (as in multipatt)
    return lev, [tuple(lev[k] for k in C[i]) for i in idx], obs, pval
