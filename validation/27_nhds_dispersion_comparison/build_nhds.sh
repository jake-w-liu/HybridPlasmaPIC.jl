#!/usr/bin/env bash
# Build + run NHDS — an EXTERNAL open-source plasma code (the New Hampshire Dispersion
# relation Solver, danielver02/NHDS, BSD-2-Clause) — to generate the kinetic whistler
# dispersion omega(k) that validation case 27 overlays our Hall-MHD oracle against.
#
# This is a genuine plasma code-to-code physics comparison (not an infrastructure
# library check): NHDS solves the full hot-plasma Vlasov-Maxwell dispersion relation
# independently of our code. The NHDS source / binary / output are GITIGNORED
# (regenerated here, never committed) so the repository stays small.
#
# Requires: git, gfortran, make (override the compiler with FC=...).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BD="$HERE/nhds_build"
mkdir -p "$BD"
cd "$BD"
if [ ! -d NHDS ]; then
    git clone --depth 1 https://github.com/danielver02/NHDS.git
fi
# strip the clone's own .git so editors don't register it as a nested repo and
# decorate its generated files (it is gitignored here and never updated).
rm -rf NHDS/.git
cd NHDS
if [ ! -x src/NHDS ]; then
    sh ./configure FC="${FC:-gfortran}"
    make
fi
# bundled beta=1 parallel whistler input -> output_whistler.in_plasma.dat
# (cols: k, theta, omega_re, omega_im(gamma), species, ...)
./src/NHDS whistler.in >/dev/null
OUT="$BD/NHDS/output_whistler.in_plasma.dat"
test -s "$OUT" && echo "OK: NHDS whistler dispersion written to $OUT"
