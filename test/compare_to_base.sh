#!/usr/bin/env bash
# Digitise the test sample with two builds of MuonCVXDDigitiser and compare the outputs hit by hit.
#
# Usage: compare_to_base.sh <base libMuonCVXDDigitiser.so> <head libMuonCVXDDigitiser.so>
#                           <sample.slcio> <work dir> <summary.md>
# The tables are appended to <summary.md>. Exit code: 0 if every configuration gives identical
# hits, 3 if any differs, anything else if digitising or comparing failed.
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
base_lib=$1 head_lib=$2 sample=$3 work=$4 summary=$5
mkdir -p "$work"

differ=0
for config in nosmear default; do
    bash "$HERE/run_digi.sh" "$base_lib" "$config" "$sample" "$work/base_$config"
    bash "$HERE/run_digi.sh" "$head_lib" "$config" "$sample" "$work/head_$config"
    rc=0
    python3 "$HERE/compare_hits.py" "$work/base_$config.csv" "$work/head_$config.csv" \
        --label "$config" --markdown "$summary" --differ-exit-code 3 || rc=$?
    case $rc in
        0) ;;
        3) differ=1 ;;
        *) echo "Comparing the $config outputs failed (exit code $rc)"; exit "$rc" ;;
    esac
done
[ "$differ" = 0 ] || exit 3
