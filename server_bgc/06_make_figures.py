#!/usr/bin/env python3
"""06_make_figures.py  Nature-journal-style figures (PNG 300dpi + PDF vector)."""
import csv, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7,
    "axes.linewidth": 0.6,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,  # editable text in Illustrator
})

R = "../results"
OUT = "../figures"
import os; os.makedirs(OUT, exist_ok=True)

def rd(path):
    with open(path) as f: return list(csv.DictReader(f, delimiter="\t"))

def save(fig, name):
    fig.savefig(f"{OUT}/{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{OUT}/{name}.pdf", bbox_inches="tight")
    plt.close(fig)

def sig_bracket(ax, x1, x2, y, p, h=0.02):
    ax.plot([x1,x1,x2,x2],[y,y+h,y+h,y],lw=0.6,c="black")
    label = "***" if p<0.001 else "**" if p<0.01 else "*" if p<0.05 else "n.s."
    ax.text((x1+x2)/2, y+h, label, ha="center", va="bottom", fontsize=7)

# ---------- Fig 1: presence/absence (single column, 89mm) ----------
s3 = {r["set"]: float(r["fraction"]) for r in rd(f"{R}/primary_ge3/bgc_summary.tsv")}
s2 = {r["set"]: float(r["fraction"]) for r in rd(f"{R}/sensitivity_ge2/bgc_summary.tsv")}
labels = ["Flagship", "Hub\n(\u22653 fam.)", "Hub\n(\u22652 fam.)"]
vals = [s3["flagship"], s3["hub_taxa"], s2["hub_taxa"]]
colors = ["#7f7f7f", "#c0392b", "#e67e22"]

fig, ax = plt.subplots(figsize=(89/25.4, 70/25.4))
ax.bar(labels, vals, color=colors, width=0.6, edgecolor="black", linewidth=0.5)
ax.set_ylabel("Fraction of genomes with\nco-located BGC")
ax.set_ylim(0, 1.15)
for i,v in enumerate(vals): ax.text(i, v+0.03, f"{v:.2f}", ha="center", fontsize=7)
sig_bracket(ax, 0, 1, 1.0, 4.75e-5)
sig_bracket(ax, 0, 2, 1.08, 1.4e-6)
save(fig, "fig1_presence_absence")

# ---------- Fig 2: density-normalized rate (single column) ----------
def rates(path, key):
    return [float(r["rate"]) for r in rd(path) if key in r["selected_via"]]
flag = rates(f"{R}/primary_ge3/density_normalized_rate.tsv", "flagship")
hub3 = rates(f"{R}/primary_ge3/density_normalized_rate.tsv", "hub")
hub2 = rates(f"{R}/sensitivity_ge2/density_normalized_rate.tsv", "hub")

fig, ax = plt.subplots(figsize=(89/25.4, 70/25.4))
bp = ax.boxplot([flag, hub3, hub2], tick_labels=labels, widths=0.5,
                patch_artist=True, showfliers=True,
                flierprops=dict(marker='o', ms=3, mfc='none', mec='black', lw=0.5),
                medianprops=dict(color="black", lw=1),
                whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6), boxprops=dict(lw=0.6))
for patch, c in zip(bp['boxes'], colors): patch.set_facecolor(c); patch.set_alpha(0.7)
ax.set_ylabel("Per-genome co-location rate")
ax.set_ylim(0, 1.25)
sig_bracket(ax, 0, 1, 1.05, 6.4e-19)
sig_bracket(ax, 0, 2, 1.15, 1.4e-6)
save(fig, "fig2_density_rate")

# ---------- Fig 3: per-genus (double column, 183mm) ----------
rows = sorted(rd(f"{R}/primary_ge3/hub_genus_breakdown.tsv"), key=lambda r:-int(r["n_genomes"]))
genera = [r["genus"].capitalize() if not r["genus"].isupper() and not r["genus"][0].isdigit()
          else r["genus"] for r in rows]
frac = [float(r["fraction_colocated"]) for r in rows]
n = [int(r["n_genomes"]) for r in rows]

fig, ax = plt.subplots(figsize=(120/25.4, 75/25.4))
ax.bar(genera, frac, color="#2980b9", width=0.65, edgecolor="black", linewidth=0.5)
ax.axhline(s3["flagship"], color="#7f7f7f", ls="--", lw=0.8, label="Flagship baseline")
for i,(f,ni) in enumerate(zip(frac,n)): ax.text(i, f+0.03, f"n={ni}", ha="center", fontsize=6)
ax.set_ylabel("Fraction with co-located BGC")
ax.set_ylim(0, 1.15)
ax.legend(frameon=False, fontsize=6, loc="upper right")
plt.xticks(rotation=30, ha="right")
save(fig, "fig3_per_genus")

print("done ->", OUT, ": fig1/2/3 as .png (300dpi) + .pdf (vector)")