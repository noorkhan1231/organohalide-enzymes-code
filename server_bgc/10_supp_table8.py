#!/usr/bin/env python3
"""10_supp_table8.py - one CSV (Supplementary Table 8 source) from a 07 results dir; prints a summary.
usage: 10_supp_table8.py NEW_RESULTS_DIR [PREVIOUS_RESULTS_DIR]"""
import csv, os, sys
R = sys.argv[1]
OLD = sys.argv[2] if len(sys.argv) > 2 else None
GE = "≥"


def rd(p):
    return list(csv.DictReader(open(p, newline=""), delimiter="\t")) if os.path.exists(p) else []


FS = {r["set"]: r for r in rd(f"{R}/final_stats.tsv")}
CM = {r["set"]: r for r in rd(f"{R}/hitcount_CMH.tsv")}
PV = {r["set"]: r for r in rd(f"{OLD}/final_stats.tsv")} if OLD else {}
EXT, ORIG = "original 4 cells + GSTDase bphK, motifgroup_SxC", "original 4 cells"
LAB = {"primary_ge3": (f"{GE}3", EXT, 10), "primary_ge3_flag4": (f"{GE}3", ORIG, 10),
       "sensitivity_ge2": (f"{GE}2", EXT, 10), "sensitivity_ge2_flag4": (f"{GE}2", ORIG, 10),
       "primary_ge3_flank5": (f"{GE}3", EXT, 5), "primary_ge3_flank20": (f"{GE}3", EXT, 20)}
rows = []
for s, (hd, cmp_, fl) in LAB.items():
    if s not in FS:
        continue
    r, base = FS[s], s.split("_flank")[0]
    v = [x["selected_via"] for x in rd(f"{R}/{base}/selected_genomes.tsv")]
    c, p = CM.get(s, {}), PV.get(s, {})
    rows.append({"section": "comparison", "set": s, "hub_definition_families": hd, "flagship_comparator": cmp_,
        "flank_kb": fl, "hub_selected": sum("hub" in x for x in v), "flagship_selected": sum("flagship" in x for x in v),
        "both_selected": sum(x == "hub+flagship" for x in v), "hub_analysed": r["hub_n"], "hub_colocated": r["hub_col"],
        "hub_fraction": r["hub_frac"], "comparator_analysed": r["fo_n"], "comparator_colocated": r["fo_col"],
        "comparator_fraction": r["fo_frac"], "fisher_OR": r["fisher_OR"], "fisher_CI_lo": r["fisher_lo"],
        "fisher_CI_hi": r["fisher_hi"], "fisher_P": r["fisher_P"], "rate_MWU_U": r["MWU_U"], "rate_MWU_P": r["MWU_P"],
        "CMH_OR": c.get("MH_common_OR", ""), "CMH_P": c.get("P", ""),
        "previous_hub_colocated": p.get("hub_col", ""), "previous_hub_analysed": p.get("hub_n", ""),
        "previous_comparator_colocated": p.get("fo_col", ""), "previous_comparator_analysed": p.get("fo_n", ""),
        "previous_fisher_OR": p.get("fisher_OR", ""), "previous_fisher_P": p.get("fisher_P", "")})
GEN = ["Streptomyces", "Pseudomonas", "Pseudomonas_E", "Burkholderia", "Mastigocoleus"]
for d in ("primary_ge3", "sensitivity_ge2"):
    T = {x["genus"].lower(): x for x in rd(f"{R}/per_genus_{d}.tsv")}
    PT = {x["genus"].lower(): x for x in rd(f"{OLD}/per_genus_{d}.tsv")} if OLD else {}
    for g in GEN:
        x, y = T.get(g.lower()), PT.get(g.lower())
        rows.append({"section": "per_genus", "set": d, "genus": g,
            "n_genomes_analysed": x["n_genomes_analysed"] if x else 0, "n_colocated": x["n_colocated"] if x else "",
            "fraction": x["fraction"] if x else "", "total_BGC_regions": x["total_BGC_regions"] if x else "",
            "n_with_BGC": x["n_with_BGC"] if x else "", "mean_rate": x["mean_rate"] if x else "",
            "median_rate": x["median_rate"] if x else "",
            "previous_n_genomes_analysed": y["n_genomes_analysed"] if y else "",
            "previous_n_colocated": y["n_colocated"] if y else "",
            "note": "" if x else "genus absent from analysed hub genomes"})
keys = []
for r in rows:
    keys += [k for k in r if k not in keys]
out = f"{R}/SuppTable8_bgc_colocation_curated9396.csv"
with open(out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=keys, lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
print("saved", out)
for r in rows:
    if r["section"] == "comparison":
        cm = "" if r["CMH_OR"] == "" else f"; CMH OR {float(r['CMH_OR']):.2f} P {float(r['CMH_P']):.3g}"
        print(f"{r['set']}: selected hub/flagship/both {r['hub_selected']}/{r['flagship_selected']}/{r['both_selected']}; "
              f"hub {r['hub_colocated']}/{r['hub_analysed']} vs {r['comparator_colocated']}/{r['comparator_analysed']}; "
              f"OR {float(r['fisher_OR']):.2f} P {float(r['fisher_P']):.3g}{cm}")
    else:
        print(f"  {r['set']} {r['genus']}: {r['n_colocated']}/{r['n_genomes_analysed']}")
