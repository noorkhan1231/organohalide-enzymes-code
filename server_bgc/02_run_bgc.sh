#!/usr/bin/env bash
# 02_run_bgc.sh — BGC detection on the selected genomes (fna_path taken from selected_genomes.tsv).
# metaSMASH if present (--streaming off, identical detection to antiSMASH), else stock antiSMASH.
# Parallel ACROSS genomes: JOBS x CPUS_PER_JOB = 8 x 4 = 32 cores (leaves ~8 for OS/other users).
set -euo pipefail

RESULTS="../results"
SELECTED="$RESULTS/selected_genomes.tsv"
OUTBASE="../metasmash_out"
JOBS=8
CPUS_PER_JOB=4
GENEFINDER="prodigal"
EXTRA_OPTS="--cb-knownclusters"    # add --minimal to go faster / detection-only

mkdir -p "$OUTBASE"

if command -v metasmash >/dev/null 2>&1; then
    TOOL="metasmash --streaming off"; echo "[info] using metaSMASH"
elif command -v antismash >/dev/null 2>&1; then
    TOOL="antismash"; echo "[info] metaSMASH not found; using antiSMASH"
else
    echo "[error] install: conda create -n metasmash -c bioconda -c conda-forge metasmash && download-antismash-databases" >&2
    exit 1
fi

run_one() {
    local line="$1" gid fna
    IFS=$'\t' read -r gid fna <<< "$line"
    [ "$fna" = "MISSING" ] && { echo "[skip] $gid no FASTA" >&2; return 0; }
    [ -f "$fna" ] || { echo "[skip] $gid FASTA not found: $fna" >&2; return 0; }
    out="$OUTBASE/$gid"
    [ -f "$out/${gid}.json" ] && { echo "[skip] $gid done" >&2; return 0; }
    mkdir -p "$out"
    # shellcheck disable=SC2086
    $TOOL --cpus "$CPUS_PER_JOB" --genefinding-tool "$GENEFINDER" \
          --output-dir "$out" --output-basename "$gid" $EXTRA_OPTS "$fna" \
          > "$out/${gid}.log" 2>&1 \
        && echo "[ok] $gid" || echo "[FAIL] $gid (see $out/${gid}.log)" >&2
}
export -f run_one; export TOOL OUTBASE CPUS_PER_JOB GENEFINDER EXTRA_OPTS

# feed NUL-separated "genome_id<TAB>fna_path" items (cols 1 and 5), skip header
tail -n +2 "$SELECTED" | awk -F'\t' '{printf "%s\t%s\0", $1, $5}' \
  | xargs -0 -P "$JOBS" -n 1 bash -c '"'"'run_one "$@"'"'"' _

echo "[done] -> $OUTBASE/<genome_id>/<genome_id>.json"
