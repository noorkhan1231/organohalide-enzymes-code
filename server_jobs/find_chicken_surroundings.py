#!/usr/bin/env python3
"""Find shotgun metagenomes of the material chickens live in (litter, bedding, farm soil, house dust) in ENA/SRA.
Run on the server (needs internet, no CPU load):   python3 find_chicken_surroundings.py [used_projects.txt]
used_projects.txt = BioProjects already in the study, one per line (optional; the four known ones are built in).
Writes: chicken_surroundings_candidates.tsv (one row per run) and chicken_surroundings_projects.tsv (one row per project)."""
import csv, io, sys, time, urllib.parse, urllib.request, statistics, collections
API = 'https://www.ebi.ac.uk/ena/portal/api/search'
FIELDS = 'run_accession,study_accession,secondary_study_accession,sample_accession,scientific_name,study_title,sample_title,isolation_source,country,collection_date,library_strategy,library_source,library_layout,instrument_platform,base_count,center_name,first_public'
USED = {'PRJEB63360', 'PRJNA1309778', 'PRJNA1308100', 'PRJNA1392098', 'PRJNA678871'}
if len(sys.argv) > 1: USED |= {l.strip() for l in open(sys.argv[1]) if l.strip()}
BASE = 'library_source="METAGENOMIC" AND library_strategy="WGS" AND instrument_platform="ILLUMINA"'
HOST = ['poultry', 'chicken', 'broiler', 'layer', 'hen', 'hatchery', 'turkey']
MAT = ['litter', 'bedding', 'farm soil', 'barn', 'house dust', 'poultry house', 'chicken house', 'coop', 'compost', 'manure', 'cafo']
Q = ['tax_tree(1670605)']                                   # taxon 'poultry litter metagenome'
for fld in ('study_title', 'sample_title', 'isolation_source'):
    for m in MAT:
        for h in HOST: Q.append('%s="*%s*" AND %s="*%s*"' % (fld, h, fld, m))
    Q.append('%s="*poultry*litter*"' % fld)
EXCL = ('faec', 'fece', 'cecum', 'caecum', 'cecal', 'caecal', 'gut', 'cloac', 'intestin', 'ileum', 'crop ', 'carcass', 'meat', 'retail', 'egg', 'rumen', 'human', 'swab')
def get(q):
    url = API + '?' + urllib.parse.urlencode({'result': 'read_run', 'query': '(%s) AND %s' % (q, BASE), 'fields': FIELDS, 'format': 'tsv', 'limit': 0})
    for k in range(4):
        try:
            return urllib.request.urlopen(url, timeout=120).read().decode()
        except Exception as e:
            time.sleep(5 * (k + 1)); err = e
    print('FAILED', q[:80], err, file=sys.stderr); return ''
rows = {}
for i, q in enumerate(Q):
    for r in csv.DictReader(io.StringIO(get(q)), delimiter='\t'): rows[r['run_accession']] = r
    time.sleep(0.4)
    if i % 20 == 0: print(i, '/', len(Q), 'queries;', len(rows), 'runs so far', file=sys.stderr)
keep = []
for r in rows.values():
    txt = ' '.join(r[k] for k in ('study_title', 'sample_title', 'isolation_source', 'scientific_name')).lower()
    if any(x in txt for x in EXCL): continue
    if r['study_accession'] in USED or r['secondary_study_accession'] in USED: continue
    try:
        if int(r['base_count']) < 1_000_000_000: continue          # under 1 Gb is too shallow
    except ValueError: continue
    keep.append(r)
keep.sort(key=lambda r: (r['study_accession'], r['run_accession']))
with open('chicken_surroundings_candidates.tsv', 'w', newline='') as o:
    w = csv.DictWriter(o, fieldnames=FIELDS.split(','), delimiter='\t', extrasaction='ignore'); w.writeheader(); w.writerows(keep)
by = collections.defaultdict(list)
for r in keep: by[r['study_accession']].append(r)
with open('chicken_surroundings_projects.tsv', 'w') as o:
    o.write('project\tsecondary\truns\tmedian_Gb\tcountry\tmaterial\ttitle\n')
    for p, v in sorted(by.items(), key=lambda kv: -len(kv[1])):
        mat = '; '.join(sorted({(x['isolation_source'] or x['sample_title'])[:40] for x in v})[:3])
        o.write('%s\t%s\t%d\t%.1f\t%s\t%s\t%s\n' % (p, v[0]['secondary_study_accession'], len(v), statistics.median(int(x['base_count']) for x in v) / 1e9,
                '; '.join(sorted({x['country'] for x in v if x['country']})[:2]), mat, v[0]['study_title'][:90]))
print('candidate runs', len(keep), 'in', len(by), 'new projects')
