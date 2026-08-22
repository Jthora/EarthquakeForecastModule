#!/usr/bin/env bash
# Cascadia tectonic tremor from the PNSN interactive tremor catalogue (Wech).
# Public, no credentials.  https://pnsn.org/tremor  |  https://tremorapi.pnsn.org
#
# Two API behaviours worth knowing:
#   - It caps a response at 20,000 events. A yearly request silently returns
#     exactly 20,000 with no error and no sign of truncation, so this fetches
#     MONTHLY and warns if any chunk approaches the cap.
#   - It returns 404, not an empty result, for a window containing no events.
#     The catalogue starts mid-2009, so early months 404 legitimately.
set -euo pipefail
OUT="$(dirname "$0")/../data/cascadia"
mkdir -p "$OUT"
CSV="$OUT/cascadia_tremor.csv"
API="https://tremorapi.pnsn.org/api/v3.0/events"
CAP=20000

START_YEAR=${1:-2009}
END_YEAR=${2:-2024}

echo "time_iso,lat,lon,depth_km,magnitude,duration_s,energy,num_stas" > "$CSV"
for y in $(seq "$START_YEAR" "$END_YEAR"); do
  for m in 01 02 03 04 05 06 07 08 09 10 11 12; do
    if [ "$m" = "12" ]; then ny=$((y+1)); nm=01; else ny=$y; nm=$(printf "%02d" $((10#$m + 1))); fi
    body=$(curl -sS --max-time 300 -w "\n%{http_code}" \
      "$API?starttime=${y}-${m}-01&endtime=${ny}-${nm}-01") || { echo "  ${y}-${m}: request failed" >&2; continue; }
    code=$(printf '%s' "$body" | tail -1)
    if [ "$code" = "404" ]; then continue; fi
    if [ "$code" != "200" ]; then echo "  ${y}-${m}: HTTP $code" >&2; continue; fi
    n=$(printf '%s' "$body" | sed '$d' \
      | python3 -c '
import json, sys, csv
w = csv.writer(sys.stderr)
d = json.load(sys.stdin)
rows = 0
for f in d.get("features", []):
    p = f.get("properties", {}); g = f.get("geometry", {}).get("coordinates", [None, None])
    if p.get("time") is None or g[0] is None:
        continue
    w.writerow([p["time"], g[1], g[0], p.get("depth"), p.get("magnitude"),
                p.get("duration"), p.get("energy"), p.get("num_stas")])
    rows += 1
print(rows)
' 2>>"$CSV")
    if [ "$n" -ge $((CAP - 1000)) ]; then
      echo "WARNING: ${y}-${m} returned $n events, near the ${CAP} cap -- likely truncated" >&2
    fi
  done
  echo "  $y done, $(( $(wc -l < "$CSV") - 1 )) rows total"
done
echo "done -> $CSV"
