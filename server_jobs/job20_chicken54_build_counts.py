"""Build FINAL-format input rows (dataset 'ambient') for the 54 new chicken surroundings runs, exactly as job1_export_counts.sh does."""
import os, glob, re, json, gzip, sys
A, O = sys.argv[1], sys.argv[2]
FAM = ['RDase','HyDase','GSTDase','OxDase','DehyDase','MetDase','FlaHase','VanHase','NHFeHase','SAMHase']
ids = [l.strip() for l in open(f'{A}/chicken54_ids.txt') if l.strip()]
miss = []
with gzip.open(f'{O}/chicken54_counts.tsv.gz', 'wt') as o, open(f'{O}/chicken54_quant_runs.tsv', 'w') as q:
    o.write('dataset\tfamily\thost_dir\trun\ttarget\tLength\tEffectiveLength\tNumReads\n')
    q.write('dataset\tfamily\thost_dir\trun\tnum_processed\tn_targets\n')
    for r in ids:
        for f in FAM:
            p = f'{A}/ambient_quant/{f}/{r}/quant.sf'
            if not os.path.exists(p): miss.append((r, f)); continue
            nproc = json.load(open(f'{os.path.dirname(p)}/aux_info/meta_info.json')).get('num_processed', '')
            ntar = 0
            with open(p) as fh:
                next(fh)
                for ln in fh:
                    c = ln.rstrip('\n').split('\t'); ntar += 1
                    if float(c[4]) > 0: o.write('\t'.join(('ambient', f, '-', r, c[0], c[1], c[2], c[4])) + '\n')
            q.write('\t'.join(map(str, ('ambient', f, '-', r, nproc, ntar))) + '\n')
def ge_of(fn):
    for l in open(fn):
        if 'genome_equivalents' in l: return re.findall(r'[0-9.eE+-]+', l.split('genome_equivalents')[-1])[-1]
    return ''
with open(f'{O}/chicken54_ge.tsv', 'w') as o:
    o.write('dataset\tfile\trun\tgenome_equivalents\n')
    for r in ids:
        fn = f'{A}/mc/{r}.txt'
        g = ge_of(fn) if os.path.exists(fn) else ''
        if not g: miss.append((r, 'GE'))
        o.write(f'ambient\t{r}\t{r}\t{g}\n')
print('runs', len(ids), 'missing', miss)
