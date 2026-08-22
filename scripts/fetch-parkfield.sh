#!/usr/bin/env bash
# Shelly's Parkfield low-frequency earthquake catalogue, from USGS ScienceBase.
# USGS-authored, public domain. 129 MB, 1,528,117 events, 88 families, 2001-2024.
set -euo pipefail
OUT="$(dirname "$0")/../data/parkfield"
mkdir -p "$OUT"
URL="https://www.sciencebase.gov/catalog/file/get/67069991d34ef5df0d802308?f=__disk__b0%2F3f%2F4f%2Fb03f4f333e893612a27c4b900311093bbaa207e0"
echo "fetching LFEcat_Apr2001-Apr2024.csv (129 MB)"
curl -sSL --fail -o "$OUT/LFEcat_Apr2001-Apr2024.csv" "$URL"
echo "done -> $OUT"
