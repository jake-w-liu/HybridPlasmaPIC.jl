#!/usr/bin/env bash
# Delegate build, provenance verification and cache publication to the Julia owner.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec "${JULIA:-julia}" --startup-file=no "$HERE/../reference_provenance.jl" hybridvpic "$@"
