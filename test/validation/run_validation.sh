#!/usr/bin/env bash
#
# Validation of the Gaudi port of MuonCVXDDigitiser against the Marlin processor it was ported from.
# See README.md in this directory.
#
# Usage: run_validation.sh [stage ...]
# Stages: geometry sim build digi dump analyse (default: all of them, in this order)
#
# Configuration through environment variables (defaults in brackets):
#   WORKDIR            work directory [./validation_work]
#   NEVENTS            number of simulated events [1000]
#   MULTIPLICITY       muons per event [10]
#   THETA_MIN/MAX      polar angle range of the muons, degrees [10/170]
#   PMIN/PMAX          momentum range of the muons, GeV [1/100]
#   SEED               simulation seed [2024]
#   MARLIN_REPO        repository of the Marlin processor [https://github.com/spg-berkeleylab/MuonCVXDDigitiser.git]
#   MARLIN_REF         commit of the Marlin processor [31bb9e2]
#   VARIANTS           port variants with fixes reverted [revert_cutondeltarays revert_all], empty to skip
#   OLD_IMAGE          image with Marlin and LCIO [ghcr.io/muoncollidersoft/mucoll-sim-ubuntu24:v2.11]
#   NEW_IMAGE          image with the Key4hep stack [ghcr.io/muoncollidersoft/mucoll-sim-ubuntu24:post_3_1_test]
#   CONTAINER_RUNTIME  docker or apptainer [docker]
#   JOBS               parallel build jobs [8]
#   ANALYSE_ARGS       extra arguments of analyse.py [--subdet VXDBarrel --layer 0 --check]
#
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
SRC=$(cd "$HERE/../.." && pwd)

WORKDIR=${WORKDIR:-$PWD/validation_work}
NEVENTS=${NEVENTS:-1000}
MULTIPLICITY=${MULTIPLICITY:-10}
THETA_MIN=${THETA_MIN:-10}
THETA_MAX=${THETA_MAX:-170}
PMIN=${PMIN:-1}
PMAX=${PMAX:-100}
SEED=${SEED:-2024}
MARLIN_REPO=${MARLIN_REPO:-https://github.com/spg-berkeleylab/MuonCVXDDigitiser.git}
MARLIN_REF=${MARLIN_REF:-31bb9e2}
VARIANTS=${VARIANTS-revert_cutondeltarays revert_all}
OLD_IMAGE=${OLD_IMAGE:-ghcr.io/muoncollidersoft/mucoll-sim-ubuntu24:v2.11}
NEW_IMAGE=${NEW_IMAGE:-ghcr.io/muoncollidersoft/mucoll-sim-ubuntu24:post_3_1_test}
CONTAINER_RUNTIME=${CONTAINER_RUNTIME:-docker}
JOBS=${JOBS:-8}
ANALYSE_ARGS=${ANALYSE_ARGS:---subdet VXDBarrel --layer 0 --check}

mkdir -p "$WORKDIR"
WORKDIR=$(cd "$WORKDIR" && pwd)

# Paths inside the containers
W=/w
V=/src/test/validation
COMPACT=$W/geometry/MAIA_v0.xml

log() { echo "[$(date +%H:%M:%S)] $*"; }

# Run a shell command in an image, with the software stack set up
in_image() {
    local image=$1
    shift
    local cmd="source /opt/setup_mucoll.sh >/dev/null 2>&1; set -euo pipefail; export PYTHONPATH=$V:\${PYTHONPATH:-}; cd $W; $*"
    case $CONTAINER_RUNTIME in
        docker)
            docker run --rm -v "$WORKDIR":$W -v "$SRC":/src:ro --entrypoint bash "$image" -c "$cmd" ;;
        apptainer)
            apptainer exec --bind "$WORKDIR":$W,"$SRC":/src:ro "docker://$image" bash -c "$cmd" ;;
        *)
            echo "Unknown CONTAINER_RUNTIME '$CONTAINER_RUNTIME'" >&2; exit 1 ;;
    esac
}

# Both implementations use the MAIA_v0 geometry of the Marlin stack: the endcap drivers changed since
stage_geometry() {
    log "Copying the MAIA_v0 geometry of $OLD_IMAGE"
    in_image "$OLD_IMAGE" "rm -rf geometry && cp -r \$K4GEO/MuColl/MAIA/compact/MAIA_v0 geometry"
    log "Exporting the tracker surfaces"
    in_image "$NEW_IMAGE" "python3 $V/export_surfaces.py $COMPACT surfaces.csv > logs_export_surfaces.log 2>&1"
}

# Detailed shower mode is needed for ddsim to store the momentum and path length of the SimTrackerHits in LCIO
stage_sim() {
    log "Simulating $NEVENTS events of $MULTIPLICITY muons"
    in_image "$OLD_IMAGE" "
        ddsim --compactFile $COMPACT --enableGun --gun.particle mu- --enableDetailedShowerMode \
              --gun.momentumMin '$PMIN*GeV' --gun.momentumMax '$PMAX*GeV' \
              --gun.thetaMin '$THETA_MIN*deg' --gun.thetaMax '$THETA_MAX*deg' --gun.distribution uniform \
              --gun.multiplicity $MULTIPLICITY --random.seed $SEED --numberOfEvents $NEVENTS \
              --outputFile sim.slcio > logs_ddsim.log 2>&1
        lcio2edm4hep sim.slcio sim.edm4hep.root > logs_lcio2edm4hep.log 2>&1"
}

stage_build() {
    log "Building the Marlin processor $MARLIN_REF"
    rm -rf "$WORKDIR/marlin_src"
    git clone --quiet "$MARLIN_REPO" "$WORKDIR/marlin_src"
    git -C "$WORKDIR/marlin_src" checkout --quiet "$MARLIN_REF"
    in_image "$OLD_IMAGE" "
        cmake -S marlin_src -B marlin_build -DCMAKE_INSTALL_PREFIX=$W/marlin_install > logs_build_marlin.log 2>&1
        cmake --build marlin_build -j$JOBS >> logs_build_marlin.log 2>&1
        cmake --install marlin_build >> logs_build_marlin.log 2>&1"

    log "Building the Gaudi port"
    in_image "$NEW_IMAGE" "
        cmake -S /src -B build_port -DBUILD_TESTING=OFF > logs_build_port.log 2>&1
        cmake --build build_port -j$JOBS >> logs_build_port.log 2>&1"

    local variant fixes
    for variant in $VARIANTS; do
        case $variant in
            revert_cutondeltarays) fixes="cutondeltarays" ;;
            revert_all) fixes="cutondeltarays pixelgrid threshold" ;;
            *) echo "Unknown variant '$variant'" >&2; exit 1 ;;
        esac
        log "Building the Gaudi port with the fixes reverted: $fixes"
        rm -rf "$WORKDIR/src_$variant"
        mkdir -p "$WORKDIR/src_$variant"
        tar -C "$SRC" --exclude ./.git --exclude './build*' --exclude ./install \
            --exclude "./$(basename "$WORKDIR")" -cf - . | tar -xf - -C "$WORKDIR/src_$variant"
        python3 "$HERE/revert_fixes.py" "$WORKDIR/src_$variant" $fixes
        in_image "$NEW_IMAGE" "
            cmake -S src_$variant -B build_$variant -DBUILD_TESTING=OFF > logs_build_$variant.log 2>&1
            cmake --build build_$variant -j$JOBS >> logs_build_$variant.log 2>&1"
    done
}

stage_digi() {
    local config variant
    for config in default nosmear; do
        log "Digitising with Marlin [$config]"
        # The Marlin processor may crash when the job terminates, after the output is complete
        in_image "$OLD_IMAGE" "
            export MARLIN_DLL=\$(echo \$MARLIN_DLL | tr : '\n' | grep -vi muoncvxddigitiser | paste -sd:):$W/marlin_install/lib/libMuonCVXDDigitiser.so
            python3 $V/make_marlin_steering.py $config $COMPACT sim.slcio digi_marlin_$config.slcio steer_marlin_$config.xml
            Marlin steer_marlin_$config.xml > logs_digi_marlin_$config.log 2>&1 || true
            grep -q 'events in 1 runs written' logs_digi_marlin_$config.log || { echo 'Marlin failed, see logs_digi_marlin_$config.log'; exit 1; }"

        for variant in port $VARIANTS; do
            log "Digitising with $variant [$config]"
            in_image "$NEW_IMAGE" "
                build_$variant/run k4run $V/runValidationDigi.py --compact $COMPACT --config $config \
                    --IOSvc.Input sim.edm4hep.root --IOSvc.Output digi_${variant}_$config.edm4hep.root \
                    > logs_digi_${variant}_$config.log 2>&1"
        done
    done
}

stage_dump() {
    local config variant
    log "Dumping the Marlin outputs"
    in_image "$OLD_IMAGE" "
        for config in default nosmear; do
            python3 $V/dump_marlin.py digi_marlin_\$config.slcio hits_marlin_\$config.csv > logs_dump_marlin_\$config.log 2>&1
        done"
    log "Dumping the Gaudi outputs"
    in_image "$NEW_IMAGE" "
        python3 $V/dump_gaudi.py digi_port_default.edm4hep.root hits_port_default.csv --sim sim_hits.csv
        python3 $V/dump_gaudi.py digi_port_nosmear.edm4hep.root hits_port_nosmear.csv
        for variant in $VARIANTS; do for config in default nosmear; do
            python3 $V/dump_gaudi.py digi_\${variant}_\$config.edm4hep.root hits_\${variant}_\$config.csv
        done; done"
}

stage_analyse() {
    log "Analysing"
    in_image "$NEW_IMAGE" "python3 $V/analyse.py $W $ANALYSE_ARGS"
}

stages=("$@")
[ ${#stages[@]} -eq 0 ] && stages=(geometry sim build digi dump analyse)
for stage in "${stages[@]}"; do
    case $stage in
        geometry|sim|build|digi|dump|analyse) "stage_$stage" ;;
        *) echo "Unknown stage '$stage'" >&2; exit 1 ;;
    esac
done
log "Done: results in $WORKDIR"
