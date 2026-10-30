#!/usr/bin/env bash
set -u
export MAMBA_ROOT_PREFIX="$HOME/.local/share/mamba"
eval "$(micromamba shell hook -s bash)" 2>/dev/null
micromamba activate antismash 2>/dev/null

gid="$1"; fna="$2"
OUTBASE="../metasmash_out"
CPUS_PER_JOB=4
GENEFINDER="prodigal"
EXTRA_OPTS="--cb-knownclusters"

if command -v metasmash >/dev/null 2>&1; then TOOL="metasmash --streaming off"
elif command -v antismash >/dev/null 2>&1; then TOOL="antismash"
else echo "[FAIL] $gid antismash not on PATH even after activation" >&2; exit 1; fi

[ "$fna" = "MISSING" ] && { echo "[skip] $gid no FASTA" >&2; exit 0; }
[ -f "$fna" ] || { echo "[skip] $gid FASTA not found: $fna" >&2; exit 0; }
out="$OUTBASE/$gid"
log="$OUTBASE/${gid}.log"
if [ -f "$out/${gid}.json" ]; then echo "[skip] $gid done" >&2; exit 0; fi
rm -rf "$out"

$TOOL --cpus "$CPUS_PER_JOB" --genefinding-tool "$GENEFINDER" \
      --output-dir "$out" --output-basename "$gid" $EXTRA_OPTS "$fna" \
      > "$log" 2>&1 \
  && echo "[ok] $gid" || echo "[FAIL] $gid (see $log)" >&2
