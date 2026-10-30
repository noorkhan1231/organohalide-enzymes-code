#!/usr/bin/env python3
"""07_final_stats_figures.py - final statistics, tables and figures for the BGC
co-location analysis. Everything is recomputed from per-genome files (no hard-coded
P values). Comparisons use disjoint groups: hub genomes vs flagship-only genomes
(a genome that is both hub and flagship counts as hub)."""
import csv, glob, os, statistics, sys
from collections import defaultdict
from scipy.stats import fisher_exact, mannwhitneyu
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

A = sys.argv[1:] + [None]*4
R, OUT = A[0] or "../results", A[1] or "../figures"
HUBCSV = A[3] or ("/home/students/yangkaiwen/Noor/Salmon_Mapping/analysis/subfamily/"
          "FINAL_RESULTS/tables/hub_genera_subfamily_summary.csv")
TAXGLOB = A[2] or "/home/students/yangkaiwen/Noor/Taxonomy_Annotation/*_taxonomy_complete.tsv"
SETS = ["primary_ge3", "sensitivity_ge2", "primary_ge3_flank5", "primary_ge3_flank20",
        "primary_ge3_flag4", "sensitivity_ge2_flag4"]
EXPECT = {"primary_ge3": (49761.5, 6.375e-19), "sensitivity_ge2": (187945.5, 1.433e-06)}
GE, X = "\u2265", "\u00d7"
C_FO, C_H3, C_H2 = "#7f7f7f", "#c0392b", "#e67e22"
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"], "font.size": 7,
    "axes.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "pdf.fonttype": 42})


def rd(p, d="\t"):
    with open(p, newline="") as f:
        return list(csv.DictReader(f, delimiter=d))


def wr(p, header, rows):
    with open(p, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def per_genome(s):
    g = {}
    for r in rd(f"{R}/{s}/bgc_per_genome.tsv"):
        x = g.setdefault(r["genome_id"], {"via": r["selected_via"], "tot": 0, "col": 0})
        if r["contig"] == "-":
            continue
        x["tot"] += 1
        if r["relationship"] == "within" or r["relationship"].startswith("flank"):
            x["col"] += 1
    return g


def split(g):
    return ({k: v for k, v in g.items() if "hub" in v["via"]},
            {k: v for k, v in g.items() if v["via"] == "flagship"})


def rates(grp):
    return [v["col"] / v["tot"] for v in grp.values() if v["tot"]]


KEYS = ["hub_n", "hub_col", "hub_frac", "fo_n", "fo_col", "fo_frac", "fisher_OR", "fisher_lo", "fisher_hi", "fisher_P",
        "hub_rate_n", "hub_rate_mean", "hub_rate_median", "fo_rate_n", "fo_rate_mean",
        "fo_rate_median", "MWU_U", "MWU_P"]


def test(g):
    hub, fo = split(g)
    hk = sum(v["col"] > 0 for v in hub.values())
    fk = sum(v["col"] > 0 for v in fo.values())
    OR, pf = fisher_exact([[hk, len(hub) - hk], [fk, len(fo) - fk]])
    import math
    se = math.sqrt(1 / hk + 1 / (len(hub) - hk) + 1 / fk + 1 / (len(fo) - fk))
    lo, hi = math.exp(math.log(OR) - 1.96 * se), math.exp(math.log(OR) + 1.96 * se)
    hr, fr = rates(hub), rates(fo)
    U, pm = mannwhitneyu(hr, fr, alternative="two-sided")
    t = dict(hub_n=len(hub), hub_col=hk, hub_frac=hk / len(hub), fo_n=len(fo), fo_col=fk,
             fo_frac=fk / len(fo), fisher_OR=OR, fisher_lo=lo, fisher_hi=hi, fisher_P=pf, hub_rate_n=len(hr),
             hub_rate_mean=statistics.mean(hr), hub_rate_median=statistics.median(hr),
             fo_rate_n=len(fr), fo_rate_mean=statistics.mean(fr),
             fo_rate_median=statistics.median(fr), MWU_U=U, MWU_P=pm)
    return t, hr, fr


# ---------------- statistics for every available set ----------------
cache, stat_rows = {}, []
for s in SETS:
    if not os.path.exists(f"{R}/{s}/bgc_per_genome.tsv"):
        print(f"[skip] {s}: {R}/{s}/bgc_per_genome.tsv not found")
        continue
    g = per_genome(s)
    t, hr, fr = test(g)
    cache[s] = (g, t, hr, fr)
    stat_rows.append([s] + [t[k] for k in KEYS])
    print(f"{s}: hub {t['hub_col']}/{t['hub_n']} vs flagship-only {t['fo_col']}/{t['fo_n']} "
          f"Fisher OR={t['fisher_OR']:.3f} P={t['fisher_P']:.4g} | rate MWU U={t['MWU_U']} "
          f"P={t['MWU_P']:.4g}")
    if s in EXPECT and "old" in R:
        ok = abs(t["MWU_U"] - EXPECT[s][0]) < 1e-6
        print(f"   reconciliation with earlier 04 output (U={EXPECT[s][0]}): "
              f"{'MATCH' if ok else 'MISMATCH - investigate before using'}")
wr(f"{R}/final_stats.tsv", ["set"] + KEYS, stat_rows)

# GSTDase-only flagship genomes (extended flagship-only minus original-4 flagship-only) and BGC totals
gst = []
for s, s4 in (("primary_ge3", "primary_ge3_flag4"), ("sensitivity_ge2", "sensitivity_ge2_flag4")):
    if s in cache and s4 in cache:
        fo, fo4 = split(cache[s][0])[1], split(cache[s4][0])[1]
        G = [v for k, v in fo.items() if k not in fo4]
        kk = sum(v["col"] > 0 for v in G)
        gst.append([s, len(G), kk, kk / len(G) if G else float("nan")])
        print(f"GSTDase-only flagship genomes ({s}): {kk}/{len(G)}")
wr(f"{R}/gstdase_only.tsv", ["set", "n_genomes", "n_colocated", "fraction"], gst)
tot = []
for s in cache:
    for nm, grp in zip(("hub", "flagship_only"), split(cache[s][0])):
        tot.append([s, nm, len(grp), sum(v["tot"] > 0 for v in grp.values()),
                    sum(v["tot"] for v in grp.values()), sum(v["col"] for v in grp.values())])
wr(f"{R}/bgc_totals.tsv", ["set", "group", "genomes", "genomes_with_region", "total_regions",
   "colocated_regions"], tot)

# ---------------- per-genus (analysed genomes only) ----------------
canon = {r["genus"].lower(): r["genus"] for r in rd(HUBCSV, ",")}


def per_genus(s):
    g = cache[s][0]
    sel = rd(f"{R}/{s}/selected_genomes.tsv")
    gm = {r["genome_id"]: r["hub_taxon"].split(";")[0] for r in sel if "hub" in r["selected_via"]}
    agg = defaultdict(list)
    unmapped = 0
    for gid, v in g.items():
        if "hub" not in v["via"]:
            continue
        if gid not in gm:
            unmapped += 1
            continue
        agg[canon.get(gm[gid], gm[gid])].append(v)
    if unmapped:
        print(f"[warn] {s}: {unmapped} hub genomes without a genus label")
    rows = []
    for gen, L in agg.items():
        rr = [v["col"] / v["tot"] for v in L if v["tot"]]
        k = sum(v["col"] > 0 for v in L)
        rows.append([gen, len(L), k, k / len(L), sum(v["tot"] for v in L), len(rr),
                     statistics.mean(rr) if rr else float("nan"),
                     statistics.median(rr) if rr else float("nan")])
    rows.sort(key=lambda r: -r[1])
    wr(f"{R}/per_genus_{s}.tsv", ["genus", "n_genomes_analysed", "n_colocated", "fraction",
       "total_BGC_regions", "n_with_BGC", "mean_rate", "median_rate"], rows)
    return rows, gm


pg = {s: per_genus(s) for s in ("primary_ge3", "sensitivity_ge2") if s in cache}

# ---------------- target-hit-count confound check ----------------
strata_rows, cmh_rows = [], []
from scipy.stats import chi2


def cmh(T):
    num = var = ar = br = 0.0
    for a, b, c, d in T:
        n = a + b + c + d
        if n < 2 or a + b == 0 or c + d == 0:
            continue
        num += a - (a + b) * (a + c) / n
        var += (a + b) * (c + d) * (a + c) * (b + d) / (n * n * (n - 1))
        ar += a * d / n
        br += b * c / n
    x = (abs(num) - 0.5) ** 2 / var if var > 0 else float("nan")
    return (ar / br if br else float("nan")), x, chi2.sf(x, 1)

for s in ("primary_ge3", "sensitivity_ge2", "primary_ge3_flag4", "sensitivity_ge2_flag4"):
    if s not in cache:
        continue
    hc = defaultdict(int)
    for r in rd(f"{R}/{s}/selected_hits_with_coords.tsv"):
        hc[r["genome_id"]] += 1
    hub, fo = split(cache[s][0])
    T = []
    for lab, lo, hi in (("0", 0, 0), ("1", 1, 1), ("2", 2, 2), (">=3", 3, 10 ** 9)):
        H = [v for k, v in hub.items() if lo <= hc.get(k, 0) <= hi]
        F = [v for k, v in fo.items() if lo <= hc.get(k, 0) <= hi]
        a = sum(v["col"] > 0 for v in H); c = sum(v["col"] > 0 for v in F)
        T.append((a, len(H) - a, c, len(F) - c))
    cmh_rows.append([s] + list(cmh(T)) + [str(T)])
    for name, grp in (("hub", hub), ("flagship_only", fo)):
        for lab, lo, hi in (("0", 0, 0), ("1", 1, 1), ("2", 2, 2), (">=3", 3, 10 ** 9)):
            L = [v for k, v in grp.items() if lo <= hc.get(k, 0) <= hi]
            k = sum(v["col"] > 0 for v in L)
            strata_rows.append([s, name, lab, len(L), k, k / len(L) if L else float("nan"), ""])
    hh, ff = [hc.get(k, 0) for k in hub], [hc.get(k, 0) for k in fo]
    U, p = mannwhitneyu(hh, ff, alternative="two-sided")
    strata_rows.append([s, "hits_per_genome", f"hub_median={statistics.median(hh)};"
                        f"fo_median={statistics.median(ff)}", len(hh), len(ff), U, p])
wr(f"{R}/hitcount_strata.tsv", ["set", "group", "resolved_hits_per_genome", "n_genomes",
   "n_colocated", "fraction_or_U", "P"], strata_rows)

wr(f"{R}/hitcount_CMH.tsv", ["set", "MH_common_OR", "chi2", "P", "strata_tables_hubcol_hubnot_focol_fonot"], cmh_rows)
for r in cmh_rows:
    print(f"hit-density-adjusted (CMH) {r[0]}: OR={r[1]:.3f} P={r[3]:.4g}")

# ---------------- genus-assignment audit (substring matching in 01) ----------------
lin = defaultdict(set)
for f in glob.glob(TAXGLOB):
    with open(f, newline="") as fh:
        rdr = csv.reader(fh, delimiter="\t")
        next(rdr, None)
        for r in rdr:
            if len(r) >= 5:
                gg = [x[3:].lower() for x in r[4].split(";") if x.startswith("g__")]
                if gg:
                    lin[r[1]].add(gg[0])
audit = []
for s, (_, gm) in pg.items():
    for gid, assigned in gm.items():
        if assigned not in lin.get(gid, set()):
            audit.append([s, gid, assigned, ";".join(sorted(lin.get(gid, set())))])
wr(f"{R}/genus_match_audit.tsv", ["set", "genome_id", "assigned_hub_genus",
   "exact_GTDB_genus"], audit)
print(f"genus audit: {len(audit)} hub genomes whose exact GTDB genus differs from the "
      f"assigned hub genus (0 expected)")


# ---------------- figures ----------------
import random


def pfmt(p):
    if p <= 0:
        return "$P$ < 10$^{-300}$"
    m, e = f"{p:.2e}".split("e")
    return rf"$P$ = {m} {X} 10$^{{{int(e)}}}$"


def bracket(ax, y, p, h):
    ax.plot([0, 0, 1, 1], [y, y + h, y + h, y], lw=0.6, c="k")
    ax.text(0.5, y + h * 1.3, pfmt(p), ha="center", va="bottom", fontsize=6)


def save(fig, name):
    fig.savefig(f"{OUT}/{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{OUT}/{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def boot(v, n=2000):
    random.seed(1)
    m = sorted(statistics.mean(random.choices(v, k=len(v))) for _ in range(n))
    return m[int(.025 * n)], m[int(.975 * n)]


def lab(name, n):
    base = name.split("_")[0]
    if base.isalpha() and base[0].isupper() and base[1:].islower():
        return rf"$\mathit{{{name.replace('_', chr(92) + '_')}}}$ (n = {n})"
    return f"{name} (n = {n})"


if "primary_ge3" in cache and "sensitivity_ge2" in cache:
    fig = plt.figure(figsize=(183 / 25.4, 125 / 25.4))
    gt = fig.add_gridspec(1, 2, left=0.08, right=0.98, top=0.96, bottom=0.58, wspace=0.35)
    gb = fig.add_gridspec(1, 1, left=0.40, right=0.97, top=0.44, bottom=0.11)
    for col, (s, s4, hl, c) in enumerate((("primary_ge3", "primary_ge3_flag4", f"Hub ({GE}3 families)", C_H3),
                                          ("sensitivity_ge2", "sensitivity_ge2_flag4", f"Hub ({GE}2 families)", C_H2))):
        t = cache[s][1]
        t4 = cache[s4][1] if s4 in cache else None
        grp = []
        if t4:
            grp.append((f"Flagship-only\n(original 4 cells)\n(n = {t4['fo_n']})", t4["fo_frac"], "#b8b8b8"))
        grp.append((f"Flagship-only\n(+GSTDase cells)\n(n = {t['fo_n']})", t["fo_frac"], C_FO))
        grp.append((f"{hl}\n(n = {t['hub_n']})", t["hub_frac"], c))
        ax = fig.add_subplot(gt[0, col])
        ax.bar(range(len(grp)), [g[1] for g in grp], color=[g[2] for g in grp], width=0.65, edgecolor="k", lw=0.5)
        top = max(g[1] for g in grp)
        for i, g in enumerate(grp):
            ax.text(i, g[1] + top * 0.02, f"{g[1]:.3f}", ha="center", fontsize=6)
        ax.set_ylim(0, top * 1.2)
        ax.set_xticks(range(len(grp)))
        ax.set_xticklabels([g[0] for g in grp], fontsize=6)
        ax.set_ylabel("Fraction of genomes with a\nco-located BGC")
        ax.text(-0.22, 1.04, "ab"[col], transform=ax.transAxes, fontsize=8, fontweight="bold")
    FR = [("primary_ge3", f"Hub {GE}3 vs flagship-only (+GSTDase cells)"),
          ("primary_ge3_flag4", f"Hub {GE}3 vs flagship-only (original 4 cells)"),
          ("sensitivity_ge2", f"Hub {GE}2 vs flagship-only (+GSTDase cells)"),
          ("sensitivity_ge2_flag4", f"Hub {GE}2 vs flagship-only (original 4 cells)"),
          ("primary_ge3_flank5", f"Hub {GE}3, 5 kb flank (+GSTDase cells)"),
          ("primary_ge3_flank20", f"Hub {GE}3, 20 kb flank (+GSTDase cells)")]
    FR = [x for x in FR if x[0] in cache][::-1]
    ax = fig.add_subplot(gb[0, 0])
    for i, (k, l) in enumerate(FR):
        t = cache[k][1]
        ax.errorbar(t["fisher_OR"], i, xerr=[[t["fisher_OR"] - t["fisher_lo"]], [t["fisher_hi"] - t["fisher_OR"]]],
                    fmt="o", ms=3, c="k", lw=0.6, capsize=2)
        ax.text(t["fisher_hi"] * 1.15, i, f"OR = {t['fisher_OR']:.2f}; {pfmt(t['fisher_P'])}", va="center", fontsize=6)
    ax.axvline(1, ls="--", lw=0.8, c="#555")
    ax.set_xscale("log")
    xmax = max(cache[k][1]["fisher_hi"] for k, _ in FR)
    ax.set_xlim(0.8, xmax * 12)
    ax.set_yticks(range(len(FR)))
    ax.set_yticklabels([l for _, l in FR][:], fontsize=6)
    ax.set_xlabel("Odds ratio, hub vs flagship-only (Fisher exact test; 95% CI, Woolf)")
    fig.text(0.01, 0.46, "c", fontsize=8, fontweight="bold")
    save(fig, "figBGC1_colocation")

if pg.get("primary_ge3", ([],))[0]:
    rows = pg["primary_ge3"][0]
    t = cache["primary_ge3"][1]
    fig, axs = plt.subplots(1, 2, figsize=(183 / 25.4, 100 / 25.4), sharey=True)
    ypos = list(range(len(rows)))[::-1]
    for ax, (vi, base, xl) in zip(axs, ((3, t["fo_frac"], "Fraction of genomes with a\nco-located BGC"),
                                       (6, t["fo_rate_mean"], "Mean co-located fraction of\nBGC regions per genome"))):
        vals = [0 if r[vi] != r[vi] else r[vi] for r in rows]
        ax.barh(ypos, vals, color=C_H3, height=0.7, edgecolor="k", lw=0.5)
        for yy, v in zip(ypos, vals):
            if v == 0:
                ax.text(max(vals + [base]) * 0.012, yy, "0", va="center", fontsize=6)
        ax.axvline(base, ls="--", lw=0.8, c="#555", label="Flagship-only (+GSTDase cells)")
        ax.set_xlim(0, max(vals + [base]) * 1.1)
        ax.set_xlabel(xl)
    h_, l_ = axs[0].get_legend_handles_labels()
    fig.legend(h_, l_, loc="lower center", frameon=False, fontsize=6, bbox_to_anchor=(0.55, -0.04))
    axs[0].set_yticks(ypos)
    axs[0].set_yticklabels([lab(r[0], r[1]) for r in rows])
    for ax, L in zip(axs, "ab"):
        ax.text(-0.05 if L == "b" else -0.55, 1.02, L, transform=ax.transAxes, fontsize=8, fontweight="bold")
    fig.tight_layout()
    save(fig, "figBGC2_per_genus")

print(f"done -> {R}/final_stats.tsv, per_genus_*.tsv, hitcount_strata.tsv, "
      f"genus_match_audit.tsv; {OUT}/figBGC1_colocation, figBGC2_per_genus (.png 300 dpi + .pdf)")
