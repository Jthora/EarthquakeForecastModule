#!/usr/bin/env bash
# Global earthquake catalogue from USGS ComCat (FDSN event service).
# Public, no credentials.  https://earthquake.usgs.gov/fdsnws/event/1/
#
# Defaults to M >= 5.5 from 1970, which the decadal completeness check found to be
# the stable threshold: counts per decade run 4377, 4384, 4865, 5160, 4898 -- 18%
# spread with no trend. M5.0+ instead climbs monotonically 13,581 -> 18,469 across
# the same decades, i.e. it is incomplete early. Mc drift projects onto long-period
# features and manufactures exactly the signal the band prediction looks for, so
# the threshold is not a detail.
#
# The service caps a response at 20,000 events, so this fetches by decade.
set -euo pipefail
OUT="$(dirname "$0")/../data/comcat"
mkdir -p "$OUT"
CSV="$OUT/global_m55.csv"
API="https://earthquake.usgs.gov/fdsnws/event/1/query"

MINMAG=${1:-5.5}
START=${2:-1970}
END=${3:-2025}

echo "time,latitude,longitude,depth,mag" > "$CSV"
for d in $(seq "$START" 10 $((END - 1))); do
  hi=$((d + 10)); [ "$hi" -gt "$END" ] && hi=$END
  echo -n "  ${d}-${hi} ... "
  curl -sS --fail --max-time 300 \
    "$API?format=csv&starttime=${d}-01-01&endtime=${hi}-01-01&minmagnitude=${MINMAG}&orderby=time-asc" \
    | python3 -c '
import sys, csv
r = csv.DictReader(sys.stdin)
w = csv.writer(sys.stdout)
n = 0
for row in r:
    try:
        w.writerow([row["time"], row["latitude"], row["longitude"],
                    row["depth"], row["mag"]])
        n += 1
    except KeyError:
        continue
print(f"{n} rows", file=sys.stderr)
' >> "$CSV"
done
echo "done -> $CSV  ($(( $(wc -l < "$CSV") - 1 )) events)"
