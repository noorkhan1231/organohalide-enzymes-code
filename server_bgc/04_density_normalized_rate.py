#!/usr/bin/env python3
"""
04_density_normalized_rate.py — controls for raw BGC density when comparing
hub-taxon vs. flagship genomes on co-location with target gene hits.

Rationale: hub_taxa (Streptomyces) genomes average far more BGCs/genome than
flagship genomes (see bgc_summary.tsv: 986/27 ~ 36.5 vs 3959/636 ~ 6.2). A
reviewer's fair objection to the presence/absence result (OR=5.48, p<0.0001)
is that hub genomes could simply have "more BGCs everywhere", inflating the
chance of ANY co-location without genes and clusters being more spatially
linked than expected by chance.

This script computes, per genome, the density-normalized rate:
    rate = (# of that genome's BGCs classified within/flanking a target hit)
           / (total BGCs detected in that genome)
and compares the hub vs. flagship distributions with Mann-Whitney U (robust
to the skew expected from BGC-count outliers), on genomes with >=1 BGC only
(rate is undefined for zero-BGC genomes).

Input:  ../results/bgc_per_genome.tsv   (from 03_bgc_gene_overlap.py)
Output: ../results/density_normalized_rate.tsv  (per-genome rate + group)
         + printed summary (median/mean per group, Mann-Whitney U, p-value)
"""
import argparse, csv, sys
from collections import defaultdict

try:
    from scipy.stats import mannwhitneyu
except ImportError:
    sys.exit("scipy not found in this env — run `pip install scipy` or use the antismash/base env that has it")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-genome", default="../results/bgc_per_genome.tsv")
    ap.add_argument("--out", default="../results/density_normalized_rate.tsv")
    args = ap.parse_args()

    # genome_id -> {via, total_bgc, colocated_bgc}
    stats = defaultdict(lambda: {"via": "-", "total": 0, "coloc": 0})

    with open(args.per_genome, newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for row in rd:
            gid = row["genome_id"]
            stats[gid]["via"] = row["selected_via"]
            if row["contig"] == "-":          # the "no_BGC_detected" placeholder row
                continue
            stats[gid]["total"] += 1
            if row["relationship"] == "within" or row["relationship"].startswith("flank"):
                stats[gid]["coloc"] += 1

    rows = []
    for gid, d in stats.items():
        if d["total"] == 0:
            continue  # rate undefined for genomes with zero BGCs
        rate = d["coloc"] / d["total"]
        rows.append((gid, d["via"], d["total"], d["coloc"], rate))

    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["genome_id", "selected_via", "total_bgc", "colocated_bgc", "rate"])
        w.writerows(sorted(rows))

    hub_rates = [r[4] for r in rows if "hub" in r[1]]
    flag_rates = [r[4] for r in rows if "hub" not in r[1] and "flagship" in r[1]]

    def summarize(name, vals):
        if not vals:
            print(f"{name}: n=0")
            return
        vals_sorted = sorted(vals)
        n = len(vals_sorted)
        median = vals_sorted[n // 2] if n % 2 else (vals_sorted[n // 2 - 1] + vals_sorted[n // 2]) / 2
        mean = sum(vals) / n
        print(f"{name}: n={n}  mean={mean:.3f}  median={median:.3f}  "
              f"min={vals_sorted[0]:.3f}  max={vals_sorted[-1]:.3f}")

    print(f"\n--- density-normalized co-location rate (BGCs with >=1 detected) ---")
    summarize("hub_taxa  (Streptomyces)", hub_rates)
    summarize("flagship  (non-hub)     ", flag_rates)

    if hub_rates and flag_rates:
        u, p = mannwhitneyu(hub_rates, flag_rates, alternative="two-sided")
        print(f"\nMann-Whitney U = {u:.1f}, p = {p:.4g}")
        if p < 0.05:
            direction = "HIGHER" if (sum(hub_rates)/len(hub_rates)) > (sum(flag_rates)/len(flag_rates)) else "LOWER"
            print(f"-> hub_taxa rate is significantly {direction} than flagship (p<0.05)")
        else:
            print("-> no significant difference in per-genome co-location RATE (p>=0.05)")
            print("   i.e. the raw presence/absence gap may be driven substantially by")
            print("   hub genomes simply carrying more BGCs overall, not tighter spatial linkage.")
    else:
        print("\nInsufficient data in one group for Mann-Whitney U test.")

    print(f"\n[done] per-genome rates written to {args.out}")


if __name__ == "__main__":
    main()
