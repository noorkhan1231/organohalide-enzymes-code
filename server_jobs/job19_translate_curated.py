#!/usr/bin/env python3
"""Clean headers and translate the 9,396 curated CDS (bacterial code 11, Prodigal start rules).
usage: python3 job19_translate_curated.py all_targets.fna curated_gene_ids.txt OUTDIR"""
import sys, re, os
fna, idf, out = sys.argv[1:4]
os.makedirs(out, exist_ok=True)
b = 'TCAG'; aa = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
T = {x+y+z: aa[i] for i, (x, y, z) in enumerate((x, y, z) for x in b for y in b for z in b)}
STARTS = {'ATG', 'GTG', 'TTG', 'CTG', 'ATT', 'ATC', 'ATA'}   # code 11 initiation codons
def read(f):
    k = h = None; d = {}
    for l in open(f):
        l = l.rstrip()
        if l.startswith('>'):
            h = l[1:]; k = h.split()[0]; d[k] = [h, []]
        elif k: d[k][1].append(l.upper())
    return {a: (v[0], ''.join(v[1])) for a, v in d.items()}
def tr(s, partial5):
    p = ''.join(T.get(s[i:i+3], 'X') for i in range(0, len(s) - 2, 3))
    if not partial5 and s[:3] in STARTS: p = 'M' + p[1:]
    return p[:-1] if p.endswith('*') else p
ids = [l.strip() for l in open(idf) if l.strip()]
seqs = read(fna)
miss = [i for i in ids if i not in seqs]
print('ids', len(ids), 'in fasta', len(seqs), 'missing', len(miss))
nstop = nM = nX = 0
with open(f'{out}/curated_genes_nucleotide.fna', 'w') as fn, open(f'{out}/curated_genes_protein.faa', 'w') as fp:
    for i in ids:
        h, s = seqs[i]
        m = re.search(r'partial=(\d)(\d)', h)
        p5 = bool(m and m.group(1) == '1')
        p = tr(s, p5)
        nstop += '*' in p; nM += p.startswith('M'); nX += 'X' in p
        for f, x in ((fn, s), (fp, p)):
            f.write('>' + i + '\n' + '\n'.join(x[j:j+80] for j in range(0, len(x), 80)) + '\n')
print('proteins with internal stop', nstop, '| start with M', nM, '| containing X (assembly N)', nX)
