#!/usr/bin/env bash
#
# Studies of the tracker digitisation in MAIAConfig: parametric (DDPlanarDigi) vs realistic
# (MuonCVXDDigitiser) hits, and the cluster shapes of the realistic digitiser.
# See README.md in this directory.
#
# Usage: run_digi_studies.sh [stage ...]
# Stages: build geometry sim digi dump analyse (default: all of them, in this order)
#
# Configuration through environment variables (defaults in brackets):
#   WORKDIR            work directory [./digi_studies_work]
#   MAIACONFIG_SRC     MAIAConfig checkout [../MAIAConfig next to this repository]
#   NEVENTS            number of simulated events [300]
#   MULTIPLICITY       muons per event [10]
#   THETA_MIN/MAX      polar angle range of the muons, degrees [0/180]
#   PMIN/PMAX          momentum range of the muons, GeV [1/100]
#   SEED               simulation seed [4242]
#   COMPACT            compact file, as seen in the image [$K4GEO/MuColl/MAIA/compact/MAIA_v0/MAIA_v0.xml]
#   IMAGE              software image [ghcr.io/muoncollidersoft/mucoll-sim-ubuntu24:post_3_1_test]
#   CONTAINER_RUNTIME  docker or apptainer [docker]
#   JOBS               parallel build jobs [8]
#
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
SRC=$(cd "$HERE/../.." && pwd)

WORKDIR=${WORKDIR:-$PWD/digi_studies_work}
MAIACONFIG_SRC=${MAIACONFIG_SRC:-$SRC/../MAIAConfig}
NEVENTS=${NEVENTS:-300}
MULTIPLICITY=${MULTIPLICITY:-10}
THETA_MIN=${THETA_MIN:-0}
THETA_MAX=${THETA_MAX:-180}
PMIN=${PMIN:-1}
PMAX=${PMAX:-100}
SEED=${SEED:-4242}
COMPACT=${COMPACT:-\$K4GEO/MuColl/MAIA/compact/MAIA_v0/MAIA_v0.xml}
IMAGE=${IMAGE:-ghcr.io/muoncollidersoft/mucoll-sim-ubuntu24:post_3_1_test}
CONTAINER_RUNTIME=${CONTAINER_RUNTIME:-docker}
JOBS=${JOBS:-8}

if [ ! -d "$MAIACONFIG_SRC/MAIAConfig/TrackerDigi" ]; then
    echo "MAIAConfig not found in '$MAIACONFIG_SRC': set MAIACONFIG_SRC to a MAIAConfig checkout" >&2
    exit 1
fi
mkdir -p "$WORKDIR"
WORKDIR=$(cd "$WORKDIR" && pwd)
MAIACONFIG_SRC=$(cd "$MAIACONFIG_SRC" && pwd)

# Paths inside the container
W=/w
D=/src/test/digi_studies
V=/src/test/validation

log() { echo "[$(date +%H:%M:%S)] $*"; }

# Run a shell command in the image, with the software stack set up
in_image() {
    local cmd="source /opt/setup_mucoll.sh >/dev/null 2>&1; set -euo pipefail; export MAIACONFIG=/maia/MAIAConfig; cd $W; $*"
    case $CONTAINER_RUNTIME in
        docker)
            docker run --rm -v "$WORKDIR":$W -v "$SRC":/src:ro -v "$MAIACONFIG_SRC":/maia:ro \
                --entrypoint bash "$IMAGE" -c "$cmd" ;;
        apptainer)
            apptainer exec --bind "$WORKDIR":$W,"$SRC":/src:ro,"$MAIACONFIG_SRC":/maia:ro \
                "docker://$IMAGE" bash -c "$cmd" ;;
        *)
            echo "Unknown CONTAINER_RUNTIME '$CONTAINER_RUNTIME'" >&2; exit 1 ;;
    esac
}

# The digitiser of this repository; build/run puts it in front of the one of the software stack
stage_build() {
    log "Building MuonCVXDDigitiser"
    in_image "cmake -S /src -B build -DBUILD_TESTING=OFF > logs_build.log 2>&1
              cmake --build build -j$JOBS >> logs_build.log 2>&1"
}

stage_geometry() {
    log "Exporting the tracker surfaces"
    in_image "PYTHONPATH=$V:\${PYTHONPATH:-} python3 $V/export_surfaces.py $COMPACT surfaces.csv > logs_surfaces.log 2>&1"
}

stage_sim() {
    log "Simulating $NEVENTS events of $MULTIPLICITY muons"
    in_image "ddsim --compactFile $COMPACT --enableGun --gun.particle mu- \
                    --gun.momentumMin '$PMIN*GeV' --gun.momentumMax '$PMAX*GeV' \
                    --gun.thetaMin '$THETA_MIN*deg' --gun.thetaMax '$THETA_MAX*deg' --gun.distribution uniform \
                    --gun.multiplicity $MULTIPLICITY --random.seed $SEED --numberOfEvents $NEVENTS \
                    --outputFile sim.edm4hep.root > logs_ddsim.log 2>&1"
}

# name and options of every digitisation
RUNS=(
    "ddplanar|"
    "ddplanar_fixed|--fix-time-window"
    "cvxd_long|--realistic --store-pixels"
    "cvxd_square|--realistic --store-pixels --tracker-pixel-v 0.025"
)

stage_digi() {
    local run name options
    for run in "${RUNS[@]}"; do
        name=${run%%|*}
        options=${run#*|}
        log "Digitising: $name"
        in_image "build/run k4run $D/tracker_digi.py $options --compact $COMPACT \
                      --IOSvc.Input sim.edm4hep.root --IOSvc.Output digi_$name.edm4hep.root > logs_digi_$name.log 2>&1"
    done
}

stage_dump() {
    log "Dumping the hits"
    in_image "python3 $D/dump_hits.py digi_ddplanar.edm4hep.root hits_ddplanar.csv sim_hits.csv
              for name in ddplanar_fixed cvxd_long cvxd_square; do
                  python3 $D/dump_hits.py digi_\$name.edm4hep.root hits_\$name.csv
              done"
}

stage_analyse() {
    log "Analysing"
    in_image "python3 -W ignore::RuntimeWarning $D/compare.py $W
              python3 -W ignore::RuntimeWarning $D/clusters.py $W"
}

stages=("$@")
[ ${#stages[@]} -eq 0 ] && stages=(build geometry sim digi dump analyse)
for stage in "${stages[@]}"; do
    case $stage in
        build|geometry|sim|digi|dump|analyse) "stage_$stage" ;;
        *) echo "Unknown stage '$stage'" >&2; exit 1 ;;
    esac
done
log "Done: results in $WORKDIR"
