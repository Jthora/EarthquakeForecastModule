#!/usr/bin/env bash
# Fetch SPICE kernels for Phase 1 (lunar tidal phase). Public, no credentials.
# de440s covers 1849-2150 at 32 MB, versus 114 MB for the full de440 — ample for
# Apollo (1969-1977) and for any terrestrial catalogue we will use.
set -euo pipefail
N="https://naif.jpl.nasa.gov/pub/naif/generic_kernels"
OUT="$(dirname "$0")/../kernels"
mkdir -p "$OUT"
fetch() { echo "fetching $(basename "$1")"; curl -sSL --fail -o "$OUT/$(basename "$1")" "$1"; }

fetch "$N/lsk/naif0012.tls"                      # leap seconds
fetch "$N/spk/planets/de440s.bsp"                # ephemeris
fetch "$N/pck/pck00011.tpc"                      # body constants
fetch "$N/pck/gm_de440.tpc"                      # GM values -- tensor is GM/d^3
fetch "$N/pck/moon_pa_de440_200625.bpc"          # lunar orientation (libration)
fetch "$N/fk/satellites/moon_de440_250416.tf"    # lunar frame definitions
echo "done -> $OUT"
