#!/usr/bin/env bash
# Run MuonCVXDDigitiser in Marlin on a sample and dump the reconstructed hits to CSV.
#
# Usage: run_digi.sh <libMuonCVXDDigitiser.so> <config> <input.slcio> <output prefix> [make_steering.py options]
#
# Only the plugins the job needs are loaded: other libraries in the stack's MARLIN_DLL crash when
# the job terminates (DD4hep teardown, MarlinKinfitProcessors), which would hide failures here.
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
lib=$1 config=$2 input=$3 out=$4
shift 4

compact=${COMPACT:-$K4GEO/MuColl/MAIA/compact/MAIA_v0/MAIA_v0.xml}
dd4hep_plugin=$(echo "${MARLIN_DLL:-}" | tr : '\n' | grep -m1 libMarlinDD4hep || true)
[ -n "$dd4hep_plugin" ] || { echo "libMarlinDD4hep not found in MARLIN_DLL"; exit 1; }
export MARLIN_DLL="$dd4hep_plugin:$lib"

python3 "$HERE/make_steering.py" "$config" "$compact" "$input" "$out.slcio" "$out.xml" "$@"
if ! Marlin "$out.xml" > "$out.log" 2>&1; then
    echo "Marlin failed, last lines of $out.log:"; tail -30 "$out.log"; exit 1
fi
grep -q "events in 1 runs written" "$out.log" || { echo "No output written, see $out.log"; exit 1; }
python3 "$HERE/dump_hits.py" "$out.slcio" "$out.csv"
