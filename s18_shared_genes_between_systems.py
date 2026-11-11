"""Step 18 (Table S5). are curated genes shared among host systems? exact-identical and >=97%-identical protein pairs between systems."""
import pandas as pd, collections, itertools, sys
from Bio import SeqIO
from Bio import Align
D=sys.argv[1].rstrip('/')+'/02_genes/'; OUT=sys.argv[2].rstrip('/')
g=pd.read_excel(D+'Supplementary_Data_3.xlsx',sheet_name='genes')
seq={r.id:str(r.seq).rstrip('*') for r in SeqIO.parse(D+'curated_genes_protein.faa','fasta')}
g['seq']=g.hit_id.map(seq); assert g.seq.notna().all(), g.seq.isna().sum()
print('genes',len(g),'genomes',g.genome_id.nunique(), 'genomes in >1 system', (g.groupby('genome_id').host_system.nunique()>1).sum())
# exact identity
ex=g.groupby('seq').agg(n=('hit_id','size'),systems=('host_system',lambda s:'|'.join(sorted(set(s)))),fam=('enzyme_family','first'))
multi=ex[ex.systems.str.contains(r'\|')]
print('distinct sequences',len(ex),'shared exactly between systems:',len(multi),'genes',multi.n.sum()); print(multi.groupby(['fam','systems']).n.agg(['size','sum']))
# >=97% identity, cross-system pairs, same family
al=Align.PairwiseAligner(); al.mode='global'; al.match_score=1; al.mismatch_score=0; al.open_gap_score=0; al.extend_gap_score=0
K=5; rows=[]
for fam,d in g.groupby('enzyme_family'):
    d=d.reset_index(drop=True); km=[set(s[i:i+K] for i in range(len(s)-K+1)) for s in d.seq]
    idx=collections.defaultdict(list)
    for i,ks in enumerate(km):
        for k in ks: idx[k].append(i)
    cnt=collections.Counter()
    for k,l in idx.items():
        if len(l)>60: continue
        for a,b in itertools.combinations(l,2):
            if d.host_system[a]!=d.host_system[b]: cnt[(a,b)]+=1
    cand=[(a,b) for (a,b),c in cnt.items() if c>=0.3*min(len(km[a]),len(km[b])) and min(len(d.seq[a]),len(d.seq[b]))>=0.9*max(len(d.seq[a]),len(d.seq[b]))]
    for a,b in cand:
        sa,sb=d.seq[a],d.seq[b]; aln=al.align(sa,sb)[0]
        ident=aln.counts().identities/max(len(sa),len(sb))
        if ident>=0.97: rows.append(dict(family=fam,gene_a=d.hit_id[a],system_a=d.host_system[a],gene_b=d.hit_id[b],system_b=d.host_system[b],identity=round(ident*100,1)))
    print(fam,'candidates',len(cand),'>=97% cross-system pairs',sum(1 for r in rows if r['family']==fam),flush=True)
R=pd.DataFrame(rows).sort_values(['family','gene_a','gene_b']).reset_index(drop=True); R.to_csv(OUT+'/cross_system_pairs_97.csv',index=False)
if len(R):
    genes=pd.unique(pd.concat([R.gene_a,R.gene_b])); print('genes in >=97% cross-system pairs',len(genes),'of',len(g),f'({100*len(genes)/len(g):.2f}%)')
    R['pair']=[ '|'.join(sorted([a,b])) for a,b in zip(R.system_a,R.system_b)]
    print(R.groupby(['family','pair']).size()); 
    gg=g.set_index('hit_id'); print(pd.Series(genes).map(gg.enzyme_family).value_counts())
