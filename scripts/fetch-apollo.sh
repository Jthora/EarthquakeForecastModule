#!/usr/bin/env bash
# Apollo Passive Seismic Experiment event catalogue, from the PDS Geosciences Node.
# Public domain, no credentials.
set -euo pipefail
B="https://pds-geosciences.wustl.edu/lunar/urn-nasa-pds-apollo_seismic_event_catalog/data"
OUT="$(dirname "$0")/../data/apollo"
mkdir -p "$OUT"
for f in levent.1008weber.csv nakamura_2005_dm_arrivals.csv nakamura_2005_dm_locations.csv; do
  echo "fetching $f"
  curl -sSL --fail -o "$OUT/$f" "$B/$f"
done
echo "done -> $OUT"
