#!/usr/bin/env bash
# Build + run Hybrid-VPIC (LANL's open-source kinetic-ion/fluid-electron hybrid PIC
# code, github.com/lanl/vpic-kokkos branch hybridVPIC) for a PERPENDICULAR shock that
# matches our setup, generating the Bz(x) profile that case 28 compares ours against.
#
# This is a live external-code, code-to-code SHOCK comparison: a different, independent
# hybrid code run from source. VPIC source/binary/output are GITIGNORED (regenerated
# here, never committed) so the repository stays small.
#
# Requires: git, cmake, and an MPI C/C++ toolchain (mpicc/mpicxx). On macOS this script
# installs OpenMPI via Homebrew if mpicxx is absent.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BD="$HERE/vpic_build"
mkdir -p "$BD"

# 1. MPI toolchain
if ! command -v mpicxx >/dev/null 2>&1; then
    if command -v brew >/dev/null 2>&1; then brew install open-mpi; else
        echo "ERROR: mpicxx not found and Homebrew unavailable; install an MPI toolchain." >&2
        exit 1
    fi
fi

# 2. clone Hybrid-VPIC (hybridVPIC branch; note: this branch does NOT use Kokkos)
cd "$BD"
[ -d vpic-kokkos ] || git clone --depth 1 -b hybridVPIC https://github.com/lanl/vpic-kokkos.git

# 3. build the VPIC deck-compiler (CMake 4.x needs the policy-min shim for old VPIC)
cd vpic-kokkos
if [ ! -x build/bin/vpic ]; then
    rm -rf build && mkdir build && cd build
    cmake -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -DCMAKE_BUILD_TYPE=Release \
          -DENABLE_INTEGRATED_TESTS=OFF -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx ..
    make -j4
    cd ..
fi

# 4. our matched perpendicular deck (committed) -> the examples dir (needs injection.cxx alongside)
cp "$HERE/shock_perp.cxx" examples/shock/

# 5. compile the deck and run it (single rank) in a clean run dir
RUN="$BD/run"; rm -rf "$RUN"; mkdir -p "$RUN"; cd "$RUN"
"$BD/vpic-kokkos/build/bin/vpic" "$BD/vpic-kokkos/examples/shock/shock_perp.cxx"
EXE=$(ls shock_perp.* | grep -v cxx | head -1)
mpirun --oversubscribe -np 1 "./$EXE"
test -s bz_profile.txt && echo "OK: Hybrid-VPIC Bz(x) profile at $RUN/bz_profile.txt"
