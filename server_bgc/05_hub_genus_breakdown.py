#!/usr/bin/env python3
"""
05_hub_genus_breakdown.py — per-genus breakdown of BGC co-location within the
hub-taxon set. Purpose: the pooled hub_taxa vs. flagship comparison treats all
hub genomes as one group; this breaks it down by genus so we can tell whether
the effect holds broadly or is carried by one or two genera.
"""
import argparse, csv, sys
from collections import defaultdict


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selected", default="../results/selected_genomes.tsv")
    ap.add_argument("--per-genome", default="../results/bgc_per_genome.tsv")
    ap.add_argument("--out", default="../results/hub_genus_breakdown.tsv")
    args = ap.parse_args()

    genome_genera = {}
    with open(args.selected, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if "hub" not in row["selected_via"]:
                continue
            genera = [g for g in row["hub_taxon"].split(";") if g and g != "-"]
            genome_genera[row["genome_id"]] = genera or ["(unmatched)"]

    gstats = defaultdict(lambda: {"total": 0, "coloc": 0})
    with open(args.per_genome, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            gid = row["genome_id"]
            if gid not in genome_genera:
                continue
            if row["contig"] == "-":
                continue
            gstats[gid]["total"] += 1
            if row["relationship"] == "within" or row["relationship"].startswith("flank"):
                gstats[gid]["coloc"] += 1

    genus_rows = defaultdict(list)
    for gid, genera in genome_genera.items():
        d = gstats.get(gid, {"total": 0, "coloc": 0})
        for g in genera:
            genus_rows[g].append((d["total"], d["coloc"]))

    out_rows = []
    for genus, vals in genus_rows.items():
        n = len(vals)
        with_bgc = [(t, c) for t, c in vals if t > 0]
        n_coloc = sum(1 for t, c in vals if c > 0)
        total_bgc = sum(t for t, c in vals)
        mean_bgc = total_bgc / n if n else 0.0
        rates = [c / t for t, c in with_bgc]
        mean_rate = sum(rates) / len(rates) if rates else float("nan")
        rates_sorted = sorted(rates)
        median_rate = (rates_sorted[len(rates_sorted)//2] if rates_sorted and len(rates_sorted) % 2
                        else (rates_sorted[len(rates_sorted)//2 - 1] + rates_sorted[len(rates_sorted)//2]) / 2
                        if rates_sorted else float("nan"))
        out_rows.append((genus, n, n_coloc, n_coloc / n if n else 0.0,
                          total_bgc, mean_bgc, mean_rate, median_rate))

    out_rows.sort(key=lambda r: -r[1])

    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["genus", "n_genomes", "n_with_colocated_BGC", "fraction_colocated",
                     "total_bgc", "mean_bgc_per_genome", "mean_rate", "median_rate"])
        for r in out_rows:
            w.writerow([r[0], r[1], r[2], f"{r[3]:.3f}", r[4], f"{r[5]:.2f}",
                        "nan" if r[6] != r[6] else f"{r[6]:.3f}",
                        "nan" if r[7] != r[7] else f"{r[7]:.3f}"])

    total_hub_n = sum(r[1] for r in out_rows)
    print(f"{'genus':<18} {'n':>4} {'w/coloc':>8} {'frac':>6} {'tot_bgc':>8} "
          f"{'mean_bgc':>9} {'mean_rate':>10} {'median_rate':>12}   %of_hub_n")
    for r in out_rows:
        pct = 100 * r[1] / total_hub_n if total_hub_n else 0
        mr = "  n/a " if r[6] != r[6] else f"{r[6]:.3f}"
        md = "  n/a  " if r[7] != r[7] else f"{r[7]:.3f}"
        print(f"{r[0]:<18} {r[1]:>4} {r[2]:>8} {r[3]:>6.3f} {r[4]:>8} "
              f"{r[5]:>9.2f} {mr:>10} {md:>12}   {pct:5.1f}%")

    strep = next((r for r in out_rows if r[0] == "Streptomyces"), None)
    if strep:
        print(f"\nStreptomyces share of pooled hub n: {strep[1]}/{total_hub_n} "
              f"({100*strep[1]/total_hub_n:.1f}%)")

    print(f"\n[done] {args.out}")


if __name__ == "__main__":
    main()
