#!/usr/bin/env python3
"""
01_select_and_resolve.py    BGC co-location

Builds the bounded genome set and resolves coordinates for every target-gene hit
inside those genomes, so 02 (antiSMASH) and 03 (overlap) need no further input.

  Set A (hub taxa):  genomes whose GTDB lineage (col 5 of *_taxonomy_complete.tsv)
                     contains a genus in hub_taxa.txt  (default: Streptomyces).
  Set B (flagship):  genomes carrying one of the four flagship gene/subfamily cells
                     in master_subfamily_assignments.csv.

For every selected genome, ALL of its target-gene hits (across the 10 gene families)
are pulled and their contig/start/end/strand resolved from the genome's GFF
(reusing the roots + regexes from your genomic_context_check resolver).

Outputs (in --outdir):
  selected_genomes.tsv         genome_id, selected_via, hub_taxon, environment, fna_path
  selected_hits_with_coords.tsv  genome_id, hit_id, gene, environment, contig, start, end, strand, selected_via
  unresolved_hits.tsv          hits whose coords/paths could not be found (for transparency)
"""
import argparse, csv, glob, os, re, sys
from collections import defaultdict

BASE = "/home/students/yangkaiwen/Noor"
TAX_GLOB = os.environ.get("TAX_TABLE") or f"{BASE}/Taxonomy_Annotation/All_families_taxonomy_final.tsv"
SUBFAM   = f"{BASE}/Salmon_Mapping/analysis/subfamily/05_master_assignments/master_subfamily_assignments.csv"

# taxonomy_complete.tsv columns: hit_id, genome_id, environment, gene, classification(lineage), method
FNA_SEARCH_ROOTS = [
    f"{BASE}/Metagenomes/Coral/MAG/RMD_New/Coral_genomes",
    f"{BASE}/Metagenomes/Sponge/RMD_genomes/output_directory/dereplicated_genomes",
    f"{BASE}/Metagenomes/Chicken/MGnify_genomes/MGnify_genomes",
    f"{BASE}/Metagenomes/Chicken/MAG/Genome",
    f"{BASE}/Metagenomes/Rhizosphere/MGnify_genomes",
    f"{BASE}/Traits_characteristics/missing_genomes",
]
GFF_CACHE_ROOTS = [
    f"{BASE}/Salmon_Mapping/gff_cache/Coral_gffs",
    f"{BASE}/Salmon_Mapping/gff_cache/Sponge_gffs",
    f"{BASE}/Salmon_Mapping/gff_cache/Chicken_gffs",
    f"{BASE}/Salmon_Mapping/gff_cache/Rhizosphere_gffs",
    f"{BASE}/Metagenomes/Coral/MAG/RMD_New/Coral_genomes",
    f"{BASE}/Metagenomes/Sponge/RMD_genomes/output_directory/dereplicated_genomes",
]
MAG_ID_RE  = re.compile(r"^(.+?)-scaffold_\d+_\d+$")
MGYG_ID_RE = re.compile(r"^(MGYG\d+)_\d+$")

# The four flagship (gene_family, subfamily_call) cells  V-BPO_ambiguous, NOT bare V-BPO
FLAGSHIP_CELLS = {
    ("RDase",   "motifgroup_CxxCxxCxxxCP"),
    ("HyDase",  "motifgroup_RxxR"),
    ("FlaHase", "cmlS"),
    ("VanHase", "V-BPO_ambiguous"),
}

_FNA, _GFF = {}, {}
for c in os.environ.get("FLAGSHIP_EXTRA","").split(","):
    if "/" in c: FLAGSHIP_CELLS.add(tuple(c.split("/")))


def build_index(roots, exts, cache):
    for root in roots:
        if not os.path.isdir(root):
            continue
        for dp, _, files in os.walk(root):
            for f in files:
                for e in exts:
                    if f.endswith(e):
                        cache.setdefault(f[:-len(e)], os.path.join(dp, f))
                        break


def genome_id_from_hit(hit):
    m = MGYG_ID_RE.match(hit)
    if m: return m.group(1)
    m = MAG_ID_RE.match(hit)
    if m: return m.group(1)
    return None


def parse_gff(path):
    """geneID -> (contig, start, end, strand). Stops at ##FASTA."""
    d = {}
    with open(path) as fh:
        for ln in fh:
            if ln.startswith("##FASTA"):
                break
            if ln.startswith("#") or "\t" not in ln:
                continue
            c = ln.rstrip("\n").split("\t")
            if len(c) < 9:
                continue
            attrs = c[8]
            gid = None
            for kv in attrs.split(";"):
                if kv.startswith("ID="):
                    gid = kv[3:]
                    break
            if gid is None:
                continue
            # keep first occurrence (gene/CDS)
            d.setdefault(gid, (c[0], int(c[3]), int(c[4]), c[6]))
    return d


def load_hub_taxa(path):
    if os.path.exists(path):
        taxa = [l.strip() for l in open(path) if l.strip() and not l.startswith("#")]
        if taxa:
            return taxa
    sys.stderr.write("[warn] hub_taxa.txt missing/empty -> using ['Streptomyces'] only\n")
    return ["Streptomyces"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="../results")
    ap.add_argument("--hub-file", default="hub_taxa.txt")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    hub = [h.lower() for h in load_hub_taxa(args.hub_file)]
    sys.stderr.write(f"[info] hub genera: {hub}\n")

    # ---- pass 1: read all taxonomy hits ----
    rows = []  # (hit_id, genome_id, env, gene, lineage)
    for tf in sorted(glob.glob(TAX_GLOB)):
        with open(tf) as fh:
            rd = csv.reader(fh, delimiter="\t"); next(rd, None)
            for r in rd:
                if len(r) < 5: continue
                rows.append((r[0], r[1], r[2], r[3], r[4]))
    if not rows:
        sys.exit(f"[error] no rows from {TAX_GLOB}")
    sys.stderr.write(f"[info] {len(rows)} total gene hits loaded\n")

    # Set A: hub-taxon genomes
    setA, hub_of = set(), {}
    for hit, g, env, gene, lin in rows:
        ll = lin.lower()
        for h in hub:
            if f";g__{h};" in ";"+ll+";":
                setA.add(g); hub_of.setdefault(g, h); break

    # Set B: flagship-subfamily genomes
    setB, cells_of = set(), defaultdict(set)
    with open(SUBFAM) as fh:
        rd = csv.DictReader(fh)
        for r in rd:
            if (r["gene_family"], r["subfamily_call"]) in FLAGSHIP_CELLS:
                g = genome_id_from_hit(r["tip_label"])
                if g:
                    setB.add(g); cells_of[g].add(f'{r["gene_family"]}/{r["subfamily_call"]}')

    selected = setA | setB
    via = {g: "+".join(v for v, s in [("hub", g in setA), ("flagship", g in setB)] if s)
           for g in selected}
    sys.stderr.write(f"[info] selected genomes: {len(selected)} "
                     f"(hub={len(setA)}, flagship={len(setB)}, both={len(setA & setB)})\n")

    # ---- indexes ----
    build_index(FNA_SEARCH_ROOTS, (".fa", ".fna", ".fasta"), _FNA)
    build_index(GFF_CACHE_ROOTS, (".gff", ".gff3"), _GFF)

    # ---- resolve coordinates for every hit in selected genomes ----
    hits_by_genome = defaultdict(list)
    for hit, g, env, gene, lin in rows:
        if g in selected:
            hits_by_genome[g].append((hit, env, gene))

    gff_cache = {}
    coord_rows, unresolved = [], []
    genome_rows = []
    for g in sorted(selected):
        fna = _FNA.get(g)
        gffp = _GFF.get(g)
        env0 = hits_by_genome[g][0][1] if hits_by_genome[g] else "-"
        genome_rows.append([g, via[g], hub_of.get(g, "-"), env0, fna or "MISSING"])
        gmap = {}
        if gffp:
            gff_cache.setdefault(gffp, None)
            if gff_cache[gffp] is None:
                try: gff_cache[gffp] = parse_gff(gffp)
                except Exception as e: gff_cache[gffp] = {}
            gmap = gff_cache[gffp]
        for hit, env, gene in hits_by_genome[g]:
            rec = gmap.get(hit)
            if rec is None and MAG_ID_RE.match(hit):        # fallback: contig = hit minus last _n
                contig = hit.rsplit("_", 1)[0]
                rec = ("__contig_only__", contig, None, None)
                # try to find any feature on that contig for the ordinal
                for gid_, v in gmap.items():
                    if v[0] == contig and gid_.endswith(hit.rsplit("_", 1)[1]):
                        rec = v; break
                else:
                    rec = None
            if rec is None:
                unresolved.append([g, hit, gene, "no_gff_match", gffp or "no_gff"])
                continue
            contig, s, e, strand = rec
            coord_rows.append([g, hit, gene, env, contig, s, e, strand, via[g]])

    # ---- write ----
    def w(path, header, data):
        with open(os.path.join(args.outdir, path), "w", newline="") as fh:
            cw = csv.writer(fh, delimiter="\t", lineterminator="\n"); cw.writerow(header); cw.writerows(data)

    w("selected_genomes.tsv",
      ["genome_id", "selected_via", "hub_taxon", "environment", "fna_path"], genome_rows)
    w("selected_hits_with_coords.tsv",
      ["genome_id", "hit_id", "gene", "environment", "contig", "start", "end", "strand", "selected_via"],
      coord_rows)
    w("unresolved_hits.tsv", ["genome_id", "hit_id", "gene", "reason", "gff"], unresolved)

    missing_fna = sum(1 for r in genome_rows if r[4] == "MISSING")
    sys.stderr.write(f"[done] {len(coord_rows)} hits resolved, {len(unresolved)} unresolved, "
                     f"{missing_fna} genomes missing FASTA -> {args.outdir}/\n")


if __name__ == "__main__":
    main()