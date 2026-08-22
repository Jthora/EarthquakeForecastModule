#!/usr/bin/env bash
# Build SPOTL (Agnew) for ocean tide loading strain. Free, 207 MB, includes
# Farrell Green's functions, global ocean tide models, and the land-sea database.
#
# Four things bite on a modern Mac, all handled here:
#   1. No Fortran compiler by default        -> brew install gcc
#   2. Homebrew gfortran cannot find libSystem -> export SDKROOT
#   3. ispand.c is K&R C, rejected by clang  -> -std=gnu89
#   4. Ocean models and the land-sea map ship as ASCII and must be converted to
#      the unformatted binaries nloadf reads. `mapcon f` takes NO filenames --
#      it reads stdin and writes lndsea.ind and lndsea.bit itself, with
#      status='new', so it fails if they already exist.
set -euo pipefail
DIR="${1:-$(dirname "$0")/../vendor/spotl}"
mkdir -p "$(dirname "$DIR")"

command -v gfortran >/dev/null || brew install gcc
export PATH=/opt/homebrew/bin:$PATH
export SDKROOT="$(xcrun --show-sdk-path)"

[ -d "$DIR" ] || {
  curl -sSL --fail -o /tmp/spotl.tar.gz "http://igppweb.ucsd.edu/~agnew/Spotl/spotl.tar.gz"
  mkdir -p "$DIR"; tar xzf /tmp/spotl.tar.gz -C "$(dirname "$DIR")"
}
cd "$DIR"

python3 - <<'PY'
s = open('src/Makefile').read()
s = s.replace("""#FTN = gfortran
#FFLAGS = -O3 -Wuninitialized -fno-f2c -fno-automatic -fno-range-check -fno-backslash""",
"""FTN = gfortran
FFLAGS = -O2 -fno-automatic -fno-range-check -fno-backslash -std=legacy -w
CC = gcc
CFLAGS = -O -c -std=gnu89 -Wno-implicit-int -Wno-deprecated-non-prototype""")
open('src/Makefile', 'w').write(s)
PY
(cd src && make)

# Land-sea database -> binary index.
[ -f lndsea.ind ] || gzcat lndsea/lndsea.ascii.gz | ./bin/mapcon f

# Ocean tide models -> binary. GOT4.7 is global; add regional models for coasts.
mkdir -p tm
for c in m2 s2 n2 k2 k1 o1 p1 q1 mf mm ssa; do
  f="tidmod/ascii/${c}.got4p7.2004.asc.gz"
  [ -f "$f" ] && [ ! -f "tm/${c}.got" ] && \
    gzcat "$f" | ./bin/modcon f "${c}.got" && mv "${c}.got" tm/ || true
done
echo "done. example:"
echo "  cd $DIR && ./bin/nloadf SITE 35.635 -120.150 500 tm/m2.got green/gr.gbaver.wef.p01.ce l"
