#!/usr/bin/env bash
# Run MuonCVXDDigitiser with k4run on a sample and dump the reconstructed hits to CSV.
#
# Usage: run_digi.sh <run> <config> <input.edm4hep.root> <output prefix> [make_options.py options]
#
# <run> is the run script of a build of this package (<build>/run), which sets up the environment
# for its Gaudi plugins. Writes <output prefix>.py (options), .log, .edm4hep.root and .csv.
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
run=$1 config=$2 input=$3 out=$4
shift 4

compact=${COMPACT:-$K4GEO/MuColl/MAIA/compact/MAIA_v0/MAIA_v0.xml}

python3 "$HERE/make_options.py" "$config" "$compact" "$input" "$out.edm4hep.root" "$out.py" "$@"
rm -f "$out.edm4hep.root"
if ! "$run" k4run "$out.py" > "$out.log" 2>&1; then
    echo "k4run failed, last lines of $out.log:"; tail -30 "$out.log"; exit 1
fi
[ -s "$out.edm4hep.root" ] || { echo "No output written, see $out.log"; exit 1; }
"$run" python3 "$HERE/dump_hits.py" "$out.edm4hep.root" "$out.csv"
