#!/usr/bin/env bash
# Ocean tide loading strain at every event site, via SPOTL nloadf.
#
# A single nloadf call is ~25 ms, so 18,310 sites x one constituent is about
# 8 minutes -- direct evaluation beats a HEALPix precompute here, and avoids
# interpolating a field that varies sharply across coastlines.
#
# Input:  a CSV of "lat,lon" (no header).
# Output: lat,lon,ee_amp,ee_pha,nn_amp,nn_pha,en_amp,en_pha  (nanostrain, degrees)
set -euo pipefail
SPOTL="${SPOTL_DIR:?set SPOTL_DIR to the built SPOTL tree}"
SITES="${1:?usage: ocean-loading-sites.sh sites.csv CONSTITUENT out.csv}"
CONST="${2:?}"
OUT="${3:?}"
export PATH=/opt/homebrew/bin:$PATH

echo "lat,lon,ee_amp,ee_pha,nn_amp,nn_pha,en_amp,en_pha" > "$OUT"
cd "$SPOTL"
n=0
while IFS=, read -r lat lon; do
  [ -z "$lat" ] && continue
  s=$(./bin/nloadf S "$lat" "$lon" 0 "tm/${CONST}.got" \
        green/gr.gbaver.wef.p01.ce l 2>/dev/null | awk '$1=="s"{print $2","$3","$4","$5","$6","$7}')
  echo "${lat},${lon},${s:-,,,,,}" >> "$OUT"
  n=$((n+1))
  [ $((n % 2000)) -eq 0 ] && echo "  $n sites" >&2
done < "$SITES"
echo "done -> $OUT ($n sites)" >&2
