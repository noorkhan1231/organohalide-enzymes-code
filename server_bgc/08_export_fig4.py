#!/usr/bin/env python3
"""08_export_fig4.py - Fig4c / Fig4d for the master manuscript. Reads results/*.tsv only."""
import csv, os, shutil, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm

RES = sys.argv[1] if len(sys.argv) > 1 else "../results"
ROOT = sys.argv[2] if len(sys.argv) > 2 else "/home/students/yangkaiwen/Noor/manuscript_figures/figures"
F4, OLD = f"{ROOT}/Fig4", f"{ROOT}/old"
for d in (F4, OLD):
    os.makedirs(d, exist_ok=True)
if "Liberation Sans" not in {f.name for f in fm.fontManager.ttflist}:
    sys.exit("Liberation Sans not found by matplotlib; install it before exporting")
for n in ("Fig4c_bgc_colocation.pdf", "Fig4d_bgc_per_genus.pdf"):
    if os.path.exists(f"{F4}/{n}"):
        shutil.move(f"{F4}/{n}", f"{OLD}/{n}")

LS = "Liberation Sans"
plt.rcParams.update({
    "font.family": LS, "font.size": 6.5, "axes.labelsize": 7, "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5, "legend.fontsize": 6.5, "axes.linewidth": 0.5,
    "xtick.major.width": 0.5, "ytick.major.width": 0.5, "xtick.major.size": 2, "ytick.major.size": 2,
    "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42,
    "mathtext.fontset": "custom", "mathtext.rm": LS, "mathtext.it": f"{LS}:italic",
    "mathtext.bf": f"{LS}:bold", "mathtext.default": "regular"})
W = 87.5 / 25.4
GE, X = "\u2265", "\u00d7"
G_LIGHT, G_MID, G_DARK = "#C8D4E3", "#9BAFC9", "#4C6A92"


def rd(p):
    with open(p, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def wcsv(p, rows, keys):
    with open(p, "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(keys)
        for r in rows:
            w.writerow([r[k] for k in keys])


def pf(p):
    if p >= 0.01:
        return rf"$\mathit{{P}}$ = {p:.2f}"
    m, e = f"{p:.1e}".split("e")
    return rf"$\mathit{{P}}$ = {m} {X} 10$^{{{int(e)}}}$"


S = {r["set"]: r for r in rd(f"{RES}/final_stats.tsv")}
f = lambda s, k: float(S[s][k])

# ---------------- Fig4c ----------------
fig = plt.figure(figsize=(W, 92 / 25.4))
ax = fig.add_axes([0.14, 0.58, 0.83, 0.27])
groups = (("primary_ge3", "primary_ge3_flag4", f"Hub {GE}3 families"),
          ("sensitivity_ge2", "sensitivity_ge2_flag4", f"Hub {GE}2 families"))
off, bw = (-0.28, 0.0, 0.28), 0.26
for g, (se, s4, lab) in enumerate(groups):
    for o, v, c in zip(off, (f(s4, "fo_frac"), f(se, "fo_frac"), f(se, "hub_frac")), (G_LIGHT, G_MID, G_DARK)):
        ax.bar(g + o, v, bw, color=c, edgecolor="k", lw=0.4)
        ax.text(g + o, v + 0.006, f"{v:.2f}", ha="center", fontsize=6.5)
ax.set_xticks([0, 1])
ax.set_xticklabels([f"{lab}\n(n = {int(float(S[se]['hub_n']))})" for se, _, lab in groups])
ax.set_ylim(0, 0.30)
ax.set_ylabel("Fraction of genomes with\nco-located BGC")
fig.legend([plt.Rectangle((0, 0), 1, 1, fc=c, ec="k", lw=0.4) for c in (G_LIGHT, G_MID, G_DARK)],
           ["Flagship-only, original 4 cells", "Flagship-only, + GSTDase cells", "Hub genomes"],
           loc="upper left", bbox_to_anchor=(0.14, 1.0), frameon=False, handlelength=1.2, borderaxespad=0)
rows = [("primary_ge3_flag4", f"{GE}3 vs original 4 cells"), ("primary_ge3", f"{GE}3 vs + GSTDase cells"),
        ("sensitivity_ge2_flag4", f"{GE}2 vs original 4 cells"), ("sensitivity_ge2", f"{GE}2 vs + GSTDase cells"),
        ("primary_ge3_flank5", f"{GE}3, 5 kb flank"), ("primary_ge3_flank20", f"{GE}3, 20 kb flank")][::-1]
ax2 = fig.add_axes([0.37, 0.13, 0.37, 0.29])
for i, (s, l) in enumerate(rows):
    o, lo, hi = f(s, "fisher_OR"), f(s, "fisher_lo"), f(s, "fisher_hi")
    ax2.plot([lo, hi], [i, i], c=G_DARK, lw=0.7)
    ax2.plot(o, i, "o", c=G_DARK, ms=2.5)
    ax2.text(1.05, i, pf(f(s, "fisher_P")), transform=ax2.get_yaxis_transform(), va="center")
ax2.axvline(1, ls="--", lw=0.5, c="k")
ax2.set_xscale("log")
ax2.set_xlim(0.8, 10)
ax2.set_xticks([1, 2, 5, 10])
ax2.set_xticklabels(["1", "2", "5", "10"])
ax2.minorticks_off()
ax2.set_ylim(-0.6, len(rows) - 0.4)
ax2.set_yticks(range(len(rows)))
ax2.set_yticklabels([l for _, l in rows])
ax2.set_xlabel("Odds ratio (95% CI)")
fig.savefig(f"{F4}/Fig4c_bgc_colocation.pdf")
plt.close(fig)
keys = ["set", "hub_n", "hub_col", "hub_frac", "fo_n", "fo_col", "fo_frac", "fisher_OR", "fisher_lo", "fisher_hi", "fisher_P"]
if not os.path.exists(f"{F4}/Fig4c_bgc_colocation_source_data.csv"):
    wcsv(f"{F4}/Fig4c_bgc_colocation_source_data.csv", [S[s] for s, _ in rows[::-1]], keys)

# ---------------- Fig4d ----------------
G = rd(f"{RES}/per_genus_primary_ge3.tsv")
base = [f("primary_ge3", "fo_frac"), f("primary_ge3", "fo_rate_mean")]


def lab(name, n):
    b = name.split("_")[0]
    if b.isalpha() and b[0].isupper() and b[1:].islower():
        return rf"$\mathit{{{name.replace('_', chr(92) + '_')}}}$ (n = {n})"
    return f"{name} (n = {n})"


fig = plt.figure(figsize=(W, 100 / 25.4))
axs = [fig.add_axes([0.38, 0.12, 0.27, 0.86]), fig.add_axes([0.71, 0.12, 0.27, 0.86])]
ypos = list(range(len(G)))[::-1]
specs = (("fraction", "Fraction of genomes\nwith co-located BGC", (0, 1.05), [0, 0.5, 1.0]),
         ("mean_rate", "Mean co-located fraction\nof BGC regions", (0, 0.084), [0, 0.04, 0.08]))
for k, (ax, (col, xl, lim, tk)) in enumerate(zip(axs, specs)):
    vals = [float(r[col]) for r in G]
    ax.barh(ypos, vals, height=0.7, color=("#4C6A92", "#9BAFC9")[k], edgecolor="k", lw=0.4)
    ax.axvline(base[k], ls="--", lw=0.5, c="k")
    for y, v in zip(ypos, vals):
        if v == 0:
            ax.text(lim[1] * 0.015, y, "0", va="center")
    ax.set_xlim(*lim)
    ax.set_xticks(tk)
    ax.set_ylim(-0.7, len(G) - 0.3)
    ax.set_xlabel(xl)
    if k == 0:
        ax.set_yticks(ypos)
        ax.set_yticklabels([lab(r["genus"], r["n_genomes_analysed"]) for r in G])
    else:
        ax.set_yticks([])
        ax.spines["left"].set_visible(False)
fig.savefig(f"{F4}/Fig4d_bgc_per_genus.pdf")
plt.close(fig)
if not os.path.exists(f"{F4}/Fig4d_bgc_per_genus_source_data.csv"):
    wcsv(f"{F4}/Fig4d_bgc_per_genus_source_data.csv", G, list(G[0].keys()))
print("done ->", F4)
