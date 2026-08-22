#!/usr/bin/env bash
# Global CMT focal mechanisms.  https://www.globalcmt.org/
# Public, no credentials. NDK format, five lines per event.
#
# The bulk file covers 1976-2020; later years come from monthly files.
set -euo pipefail
OUT="$(dirname "$0")/../data/gcmt"
mkdir -p "$OUT"
B="https://www.ldeo.columbia.edu/~gcmt/projects/CMT/catalog"
NDK="$OUT/gcmt.ndk"

echo "fetching jan76_dec20.ndk (23 MB)"
curl -sSL --fail --max-time 600 "$B/jan76_dec20.ndk" -o "$NDK"

MONTHS=(jan feb mar apr may jun jul aug sep oct nov dec)
for y in 2021 2022 2023 2024; do
  yy=${y:2:2}
  for m in "${MONTHS[@]}"; do
    if curl -sSL --fail --max-time 120 "$B/NEW_MONTHLY/$y/${m}${yy}.ndk" >> "$NDK" 2>/dev/null; then
      echo -n "."
    fi
  done
  echo " $y"
done
echo "done -> $NDK  ($(( $(wc -l < "$NDK") / 5 )) events)"
