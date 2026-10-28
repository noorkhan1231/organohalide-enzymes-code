# Analysis code: halogenase and dehalogenase enzymes in host-associated metagenomes

Everything for the manuscript is produced from `inputs/` by the scripts in `scripts/`, in the order below. Results go to `results/`, figures to `figures/`, tables and Supplementary Data to `tables/`, logs to `logs/`. `RESULTS_FINAL.md` holds every number; `MASTER_DATA_LEDGER.md` holds the data decisions; `ANALYSIS_STEPS.md` the step list.

## Order (run from this folder)
| Step | Script | Output |
|---|---|---|
| 1 | `scripts/s01_integrity.py inputs results/s01_integrity` | input integrity |
| 2 | `scripts/s02_metadata.py` | run metadata (ENA, reads, genome equivalents, fastp) |
| 3 | `scripts/s03_qc.py` (+ `qc_decisions.tsv`, `study_groups.tsv`), `s03b_tables.py` | samples.tsv, run_qc.tsv |
| 4 | `scripts/s04_build.py inputs results/s03_qc results/s04_build` | RPKG tables (family, subfamily, genus) |
| 5 | `scripts/s05_family.py` (+ `vegan_tests.R`, `lib/multipatt.py`) | family tests, main and GE >= 20/50/200 |
| 6 | `scripts/s06_subfamily.py`, `scripts/s06b_subfamily_nmds.R results/s06_subfamily` | subfamily tests, indicators, lineage control, NMDS |
| 6c | `Rscript scripts/s06c_subfamily_combined.R <this folder>` | combined subfamily PERMANOVA, PERMDISP, PCoA (Figure 4A) |
| 7 | `scripts/s07_metatx.py inputs results/s07_metatx` (+ `metatx_qc_decisions.tsv`) | metatranscriptomes |
| 8 | `scripts/s08_coordinates.py`, `scripts/s09_map.py` | sample coordinates, map |
| 9 | `Rscript scripts/s10_geography.R` | Mantel, latitude |
| 10 | `scripts/s11_ambient.py inputs results/s03_qc results/s04_build results/s11_ambient` | host vs surroundings with null genes |
| 11 | server: `scripts/server_bgc/` (antiSMASH 6.1.1; hub lists `hub3_genomic.txt`, `hub2_genomic.txt`) | `inputs/bgc_genomic_hubs` → `results/s12_bgc` |
| 12 | `scripts/structures/build_dist.py`, `perm.R`; `scripts/s13_structures.py` | 529 structures, structural PERMANOVA, pockets |
| 13 | `scripts/s14_hitlevel.py inputs/hitlevel inputs/bundle results/s14_hitlevel` | genes, prevalence, identity, domains, topology, plasmid contigs |
| 14 | `scripts/figures/fig1_atlas.py` … `fig6_hubs_bgc.py`, `supp_figures.py` (argument: this folder; supp also the 21-Oct package for Supp Figs 1–3) | Figures 1–6, Supplementary Figures 1–14 |
| 15 | `scripts/s15_tables.py <this folder> <21-Oct package> tables` | Table 1, Supplementary Tables 1–16, Supplementary Data 1–7 |
| 16 | `scripts/s99_verify.py <this folder>` | consistency check (all PASS on 28 Oct) |

Shared code: `scripts/stats_utils.py` (tests, BH FDR, Dunn letters), `scripts/figures/figstyle.py` (fonts, host colours, panel assembly; Liberation Sans, Chicken #E69F00, Coral #56B4E9, Rhizosphere #009E73, Sponge #CC79A7; 174 mm width).

## Software used here
Python 3.13 (pandas 3.0.5, NumPy 2.5.3, SciPy 1.18.1, statsmodels 0.15.0, matplotlib 3.11.2, GeoPandas 0.14.4, PyMuPDF 1.28.2), R 4.3.3 with vegan 2.6.4. Server tools and versions: Supplementary Table 16.
