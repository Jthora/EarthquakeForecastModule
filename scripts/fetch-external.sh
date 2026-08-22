#!/usr/bin/env bash
# Tier 2 external time series. All free, no credentials, ~200 MB total.
#
# OMNI2 is the important one: a single hourly file carrying solar wind, the
# interplanetary magnetic field and the geomagnetic indices (Kp, Dst, AE)
# together, back to 1963.
set -euo pipefail
OUT="$(dirname "$0")/../data/external"
mkdir -p "$OUT"

fetch() {
  local url="$1" name="$2" min="$3"
  echo -n "  $name ... "
  curl -sSL --fail --max-time 900 -o "$OUT/$name" "$url"
  local n; n=$(wc -c < "$OUT/$name")
  if [ "$n" -lt "$min" ]; then
    echo "SUSPECT: $n bytes, expected at least $min" >&2
  else
    echo "$((n / 1024)) KB"
  fi
}

fetch "https://spdf.gsfc.nasa.gov/pub/data/omni/low_res_omni/omni2_all_years.dat" \
      "omni2_all_years.dat" 150000000
fetch "https://kp.gfz-potsdam.de/app/files/Kp_ap_Ap_SN_F107_since_1932.txt" \
      "kp_ap_sn_f107.txt" 4000000
fetch "https://datacenter.iers.org/data/csv/finals2000A.all.csv" \
      "iers_finals2000A.csv" 3000000
fetch "https://www.sidc.be/SILSO/DATA/SN_d_tot_V2.0.csv" \
      "silso_sunspot_daily.csv" 2000000
echo "done -> $OUT"
