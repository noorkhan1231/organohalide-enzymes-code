#!/usr/bin/env bash
# Rerun BGC co-location (01->03->07->10) with the curated hit table; previous results untouched.
set -euo pipefail
cd "$(dirname "$0")"
export MAMBA_ROOT_PREFIX="$HOME/.local/share/mamba"
eval "$(micromamba shell hook -s bash)"; micromamba activate antismash
export TAX_TABLE=/home/students/yangkaiwen/Noor/Taxonomy_Annotation/All_families_taxonomy_curated9396.tsv
OUT=../results_curated9396
[ -s "$TAX_TABLE" ] || { echo "missing $TAX_TABLE"; exit 1; }
head -1 "$TAX_TABLE" | tr -d '\r' | grep -qx $'hit_id\tgenome_id\tenvironment\tgene\tclassification\tmethod' || { echo "unexpected header:"; head -1 "$TAX_TABLE"; exit 1; }
echo "== environment labels"; cut -f3 "$TAX_TABLE" | tail -n +2 | sort | uniq -c
echo "== hits per gene"; cut -f4 "$TAX_TABLE" | tail -n +2 | sort | uniq -c
mkdir -p "$OUT"
python3 - <<'PY'
s = open("01_select_and_resolve.py").read()
a = 'TAX_GLOB = f"{BASE}/Taxonomy_Annotation/All_families_taxonomy_final.tsv"'
if "TAX_TABLE" not in s:
    assert a in s, "01 TAX_GLOB line not found"
    open("01_select_and_resolve.py", "w").write(s.replace(a, 'TAX_GLOB = os.environ.get("TAX_TABLE") or ' + a[len("TAX_GLOB = "):]))
PY
export FLAGSHIP_EXTRA="GSTDase/bphK,GSTDase/motifgroup_SxC"
python 01_select_and_resolve.py --hub-file hub2.txt --outdir $OUT/sensitivity_ge2 2>&1 | tee $OUT/01_sensitivity_ge2.log
python 01_select_and_resolve.py --hub-file hub3.txt --outdir $OUT/primary_ge3 2>&1 | tee $OUT/01_primary_ge3.log
unset FLAGSHIP_EXTRA
python 01_select_and_resolve.py --hub-file hub3.txt --outdir $OUT/primary_ge3_flag4 2>&1 | tee $OUT/01_primary_ge3_flag4.log
python 01_select_and_resolve.py --hub-file hub2.txt --outdir $OUT/sensitivity_ge2_flag4 2>&1 | tee $OUT/01_sensitivity_ge2_flag4.log
cat $OUT/*/selected_genomes.tsv | grep -v '^genome_id' | cut -f1,5 | sort -u > $OUT/all_selected.tsv
awk -F'\t' '$2!="MISSING"{print $1, $2}' $OUT/all_selected.tsv | xargs -P 8 -n 2 ./run_one.sh > $OUT/run_bgc.log 2>&1 || true
echo "antiSMASH: ok $(grep -c '^\[ok\]' $OUT/run_bgc.log) skip $(grep -c '^\[skip\]' $OUT/run_bgc.log) fail $(grep -c '^\[FAIL\]' $OUT/run_bgc.log)"
nm=$(awk -F'\t' '$2!="MISSING"{print $1}' $OUT/all_selected.tsv | while read g; do [ -f ../metasmash_out/$g/$g.json ] || echo $g; done | wc -l)
[ "$nm" -eq 0 ] || { echo "$nm genomes without antiSMASH JSON (see $OUT/run_bgc.log); stopping"; exit 1; }
for s in primary_ge3 sensitivity_ge2 primary_ge3_flag4 sensitivity_ge2_flag4; do
  python 03_bgc_gene_overlap.py --selected $OUT/$s/selected_genomes.tsv --coords $OUT/$s/selected_hits_with_coords.tsv --outdir $OUT/$s
done
for k in 5 20; do
  python 03_bgc_gene_overlap.py --selected $OUT/primary_ge3/selected_genomes.tsv --coords $OUT/primary_ge3/selected_hits_with_coords.tsv --flank-kb $k --outdir $OUT/primary_ge3_flank$k
done
python 07_final_stats_figures.py $OUT ../figures_curated9396 "$TAX_TABLE"
python 10_supp_table8.py $OUT ../results
