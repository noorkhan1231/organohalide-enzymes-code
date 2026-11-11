# Halogenase and dehalogenase enzymes in host-associated metagenomes: analysis code (version 1.0.2, 10 Nov 2026)

Everything for the manuscript is produced from `inputs/` by the scripts in `scripts/`, in the order below. Results go to `results/`, figures to `figures/`, tables and Supplementary Data to `tables/`, logs to `logs/`. In this repository the scripts are in the top-level folder and in `figures/`, `server_bgc/`, `server_jobs/`, `structures/` and `lib/`; the table below calls them `scripts/...`.

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
| 10b | `scripts/chicken_lopo.py <this folder>` | chicken gut vs surroundings, each surroundings project left out in turn (Table S10E) |
| 10c | server: `server_jobs/find_chicken_surroundings.py`, `ncbi_two_stage.sh`, `process_as_ready.sh`, `null_sysmatched_chicken54.sh`, `job20_chicken54_build_counts.py`, then `setup_s11_chicken54.sh` and a rerun of step 10 | the 54 further chicken surroundings runs (round 2); the count files are in the data archive, 05_analysis_results/s11_ambient |
| 11 | server: `scripts/server_bgc/` (antiSMASH 6.1.1; hub lists `hub3_genomic.txt`, `hub2_genomic.txt`) | `inputs/bgc_genomic_hubs` → `results/s12_bgc` |
| 12 | `scripts/structures/build_dist.py`, `perm.R`; `scripts/s13_structures.py` | 529 structures, structural PERMANOVA, pockets |
| 13 | `scripts/s14_hitlevel.py inputs/hitlevel inputs/bundle results/s14_hitlevel` | genes, prevalence, identity, domains, topology, plasmid contigs |
| 13b | `scripts/s16_indicator_host_support.py inputs results results/s16_indicator_support` | indicator subfamilies vs surroundings and in metatranscriptomes (Figure 4E, Table S9) |
| 14 | `scripts/figures/fig1_atlas.py` … `fig6_hubs_bgc.py`, `supp_figures.py`, `supp_fig11_lineages.py` (argument: this folder; supp also the 21-Oct package for Supp Figs 1–3) | Figures 1–6, Supplementary Figures 1–16 (`supp_fig11_lineages.py` makes Supplementary Figure 11; the former Figures 11–15 are now 12–16) |
| 15 | `scripts/s15_tables.py <this folder> <21-Oct package> tables` | Table 1, Supplementary Tables 1–16, Supplementary Data 1–7 |
| 15b | `scripts/supp_tables_round2.py <workbook from step 15> <workbook out> <data archive folder>` | adds Table S4 B–C, S5 B–C and S10 E after steps 17–19 |
| 15c | `scripts/supp_data_add_column_definitions.py <folder with Supplementary Data workbooks> <output folder>` | adds the Column_definitions sheet (meaning of every column and abbreviation) to Supplementary Data 1, 3, 4, 5, 6, 7 |
| 16 | `scripts/s99_verify.py <this folder>` | consistency check |
| 17 | `scripts/s17_lineage_broader_ranks.py <this folder> results/s17_lineage_broader_ranks` | lineage control at GTDB family, class and phylum (Table S12) |
| 18 | `scripts/s18_shared_genes_between_systems.py <data archive folder> <output folder>` | genes shared between host systems (Table S5; about 7 min on one core) |
| 19 | `scripts/s19_reference_free_clades.py <data archive folder> <output folder>` | positions of unlabeled genes in the gene trees, reference-free clades (Table S4) |

Shared code: `scripts/stats_utils.py` (tests, BH FDR, Dunn letters), `scripts/figures/figstyle.py` (fonts, host colors, panel assembly; Liberation Sans, Chicken #E69F00, Coral #56B4E9, Rhizosphere #009E73, Sponge #CC79A7; 174 mm width).

## Software used here
Python 3.13 (pandas 3.0.5, NumPy 2.5.3, SciPy 1.18.1, statsmodels 0.15.0, matplotlib 3.11.2, GeoPandas 0.14.4, PyMuPDF 1.28.2), R 4.3.3 with vegan 2.6.4. Server tools and versions: Supplementary Table 16.

Numbering note: supplementary Tables, Figures and Data are numbered by first citation in the manuscript. Scripts keep their analysis-order internal keys and write the files under the published numbers.


`inputs/halocycdb_core_references.csv` (read by step 15) lists the 217 characterized reference enzymes shown in Table S1.

## Version history

- 1.0.2 (10 Nov 2026): chicken surroundings extended to 111 runs from 10 projects (steps 10b, 10c); steps 17–19 and `supp_tables_round2.py` added; Supplementary Figure 11 (lineage cladogram; former Figures 11–15 renumbered 12–16); scripts and comments use US spelling; absolute paths removed from `s15_tables.py` and the new scripts.
- 1.0.1 (30 Oct 2026): numbering of supplementary items by first citation.
- 1.0.0 (28 Oct 2026): first release.
