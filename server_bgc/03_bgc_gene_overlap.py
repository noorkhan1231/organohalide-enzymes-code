#!/usr/bin/env python3
"""
03_bgc_gene_overlap.py — does a detected BGC overlap / flank a target gene hit?

Consumes:
  ../metasmash_out/<gid>/<gid>.json            (from 02)
  ../results/selected_hits_with_coords.tsv     (from 01)
Emits (to ../results):
  bgc_per_genome.tsv   one row per detected BGC + nearest target hit + relationship + distance
  bgc_summary.tsv      fraction of hub / flagship / all-selected genomes with >=1 co-located BGC
Relationship: within (overlap) | flank_<=Nkb | distant. Contig must match between the
FASTA header antiSMASH read and the 'contig' column resolved from the GFF.
"""
import argparse, csv, glob, json, os, re, sys
from collections import defaultdict

LOC_RE = re.compile(r"\[?<?(\d+)>?:<?(\d+)>?\]?")


def parse_regions(jp):
    with open(jp) as fh: data = json.load(fh)
    out = []
    for rec in data.get("records", []):
        desc = (rec.get("description") or "").strip()
        rid = rec.get("id") or rec.get("name") or ""
        renamed = rid.endswith("..") or (rid.startswith("c") and "_" in rid)
        contig = desc if (desc and renamed) else (rid or desc)
        got = False
        for f in rec.get("features", []):
            if f.get("type") != "region": continue
            m = LOC_RE.search(str(f.get("location", "")))
            if not m: continue
            prods = f.get("qualifiers", {}).get("product", [])
            out.append((contig, int(m.group(1)), int(m.group(2)), ",".join(prods) or "unknown"))
            got = True
        if not got:
            for a in rec.get("areas", []):
                if a.get("start") is None: continue
                p = a.get("products") or []
                out.append((contig, int(a["start"]), int(a["end"]),
                            ",".join(p) if isinstance(p, list) else str(p)))
    return out


def relate(bgc, hits, flank):
    c, bs, be = bgc[0], bgc[1], bgc[2]
    best = ("distant", None, None); bd = None
    for hc, hs, he, hl in hits:
        if hc != c: continue
        if hs < be and he > bs: return ("within", 0, hl)
        d = hs - be if hs >= be else bs - he
        if bd is None or d < bd:
            bd = d
            best = (f"flank_<={flank//1000}kb" if d <= flank else "distant", d, hl)
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json-base", default="../metasmash_out")
    ap.add_argument("--coords", default="../results/selected_hits_with_coords.tsv")
    ap.add_argument("--selected", default="../results/selected_genomes.tsv")
    ap.add_argument("--outdir", default="../results")
    ap.add_argument("--flank-kb", type=float, default=10.0)
    args = ap.parse_args()
    flank = int(args.flank_kb * 1000)

    coords = defaultdict(list)
    with open(args.coords) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            coords[r["genome_id"]].append(
                (r["contig"], int(r["start"]), int(r["end"]),
                 f'{r["gene"]}:{r["hit_id"]}'))

    via = {}
    with open(args.selected) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            via[r["genome_id"]] = r["selected_via"]

    per, coloc, nbgc = [], {}, defaultdict(int)
    for jp in sorted(glob.glob(os.path.join(args.json_base, "*", "*.json"))):
        gid = os.path.basename(os.path.dirname(jp))
        if gid not in via:
            continue  # stale JSON from a different/broader selection run
        regs = parse_regions(jp); nbgc[gid] = len(regs)
        hits = coords.get(gid, []); c = False
        for b in regs:
            rel, d, hl = relate(b, hits, flank)
            if rel == "within" or rel.startswith("flank"): c = True
            per.append([gid, via.get(gid, "-"), b[0], b[1], b[2], b[3],
                        hl or "-", rel, "" if d is None else d])
        if not regs:
            per.append([gid, via.get(gid, "-"), "-", "-", "-", "no_BGC", "-", "-", ""])
        coloc[gid] = c

    os.makedirs(args.outdir, exist_ok=True)
    with open(os.path.join(args.outdir, "bgc_per_genome.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["genome_id", "selected_via", "contig", "bgc_start", "bgc_end",
                    "product_class", "nearest_target_hit", "relationship", "distance_bp"])
        w.writerows(per)

    sets = {"all_selected": lambda v: True,
            "hub_taxa":     lambda v: "hub" in v,
            "flagship":     lambda v: "flagship" in v}
    with open(os.path.join(args.outdir, "bgc_summary.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["set", "n_genomes", "n_with_colocated_BGC", "fraction", "total_BGCs", "flank_kb"])
        for name, fn in sets.items():
            gids = [g for g in coloc if fn(via.get(g, ""))]
            n = len(gids); k = sum(coloc[g] for g in gids)
            tot = sum(nbgc[g] for g in gids)
            w.writerow([name, n, k, f"{k/n:.3f}" if n else "nan", tot, args.flank_kb])

    sys.stderr.write(f"[done] {args.outdir}/bgc_per_genome.tsv + bgc_summary.tsv\n")


if __name__ == "__main__":
    main()
