"""Step 19 (Table S5). where do unlabeled curated genes sit in the gene trees? reference-free clades and distance to the nearest reference/outgroup tip."""
from Bio import Phylo
import pandas as pd, numpy as np, glob, re, sys
D=sys.argv[1].rstrip('/')+'/'; OUT=sys.argv[2].rstrip('/')
g=pd.read_excel(D+'02_genes/Supplementary_Data_3.xlsx',sheet_name='genes').set_index('hit_id')
hitlike=re.compile(r'^(MGYG\d+_\d+|.*_MAG_.*|.*-scaffold_.*|.*_bin.*|.*SAM[NE]A?\d+.*)$')
MIN=10; rows=[]; gene_rows=[]; clade_rows=[]
for f in sorted(glob.glob(D+'03_trees/Supplementary_Data_2/*.nwk')):
    fam=f.split('/')[-1].replace('_tree.nwk','')
    t=Phylo.read(f,'newick')
    # classify tips
    for c in t.get_terminals():
        n=c.name; c.kind='cur' if n in g.index else ('hit' if hitlike.match(n) else 'ref')
    # postorder counts
    for cl in t.find_clades(order='postorder'):
        if cl.is_terminal(): cl.nref=int(cl.kind=='ref'); cl.ncur=int(cl.kind=='cur'); cl.down=0.0 if cl.kind=='ref' else np.inf
        else:
            cl.nref=sum(k.nref for k in cl.clades); cl.ncur=sum(k.ncur for k in cl.clades)
            cl.down=min(k.down+(k.branch_length or 0) for k in cl.clades)
    # preorder: nearest ref distance from any node = min(down, up)
    t.root.up=np.inf
    for cl in t.find_clades(order='preorder'):
        for k in cl.clades:
            others=[o.down+(o.branch_length or 0) for o in cl.clades if o is not k]
            k.up=(k.branch_length or 0)+min([cl.up]+others)
    # maximal reference-free clades with >= MIN curated genes
    infree=set(); cid=0
    def walk(cl):
        global cid
        if cl.nref==0 and cl.ncur>=MIN:
            cid+=1; tips=[x.name for x in cl.get_terminals() if x.kind=='cur']
            sysc=g.loc[tips,'host_system'].value_counts().to_dict(); lab=(g.loc[tips,'subfamily']!='unassigned').sum()
            clade_rows.append(dict(family=fam,clade=f'{fam}_{cid}',curated_genes=len(tips),labeled=int(lab),systems='|'.join(sorted(sysc)),**{f'n_{k}':v for k,v in sysc.items()}))
            infree.update(tips); return
        for k in cl.clades: walk(k)
    walk(t.root)
    for c in t.get_terminals():
        if c.kind!='cur': continue
        sf=g.at[c.name,'subfamily']
        gene_rows.append(dict(family=fam,hit_id=c.name,unlabeled=sf=='unassigned' or str(sf).startswith('no_'),in_ref_free_clade=c.name in infree,dist_nearest_ref=min(c.down,c.up),identity=g.at[c.name,'identity_to_nearest_HaloCycDB_reference_pct']))
G=pd.DataFrame(gene_rows); C=pd.DataFrame(clade_rows).fillna(0)
print(g.subfamily.value_counts().head(8))
print('curated genes in trees',len(G),'unlabeled',G.unlabeled.sum())
s=G.groupby(['family','unlabeled']).agg(genes=('hit_id','size'),in_ref_free=('in_ref_free_clade','sum'),median_dist=('dist_nearest_ref','median'),median_identity=('identity','median')).reset_index()
s['pct_ref_free']=(100*s.in_ref_free/s.genes).round(1); print(s.to_string(index=False))
tot=G.groupby('unlabeled').agg(genes=('hit_id','size'),in_ref_free=('in_ref_free_clade','sum'),median_dist=('dist_nearest_ref','median'),median_identity=('identity','median')); tot['pct']=100*tot.in_ref_free/tot.genes; print(tot)
print('ref-free clades >=10 curated genes:',len(C),'; genes in them',C.curated_genes.sum(),'; single-system clades',(~C.systems.str.contains(r'\|')).sum())
print(C.groupby('family').agg(clades=('clade','size'),genes=('curated_genes','sum'),largest=('curated_genes','max'),single_system=('systems',lambda s:(~s.str.contains(r'\|')).sum())))
G.to_csv(OUT+'/gene_positions.csv',index=False); C.to_csv(OUT+'/reference_free_clades.csv',index=False); s.to_csv(OUT+'/summary_by_family.csv',index=False)
