# Validation against the Marlin processor

This directory compares the Gaudi port of `MuonCVXDDigitiser` with the Marlin processor it was ported
from ([spg-berkeleylab/MuonCVXDDigitiser](https://github.com/spg-berkeleylab/MuonCVXDDigitiser) at
31bb9e2), on the same simulated events. It is not part of `ctest`: it needs two software stacks, one
with Marlin and LCIO and one with Key4hep, and takes about half an hour.

## What is compared

A muon-gun sample is simulated with `ddsim` in the Marlin stack and converted to EDM4hep. All the MAIA
tracker sub-detectors are then digitised, with the fired pixels stored, by:

| Name | Implementation |
|---|---|
| `marlin` | The Marlin processor |
| `port` | The Gaudi port, as in this repository |
| `revert_cutondeltarays` | The Gaudi port with the `CutOnDeltaRays` fix reverted |
| `revert_all` | The Gaudi port with all the fixes that change the results reverted (`CutOnDeltaRays`, tracker barrel ladder length, threshold smearing) |

The fixes are reverted in a copy of the sources by `revert_fixes.py`. `revert_all` must agree with
Marlin: this shows that the port is faithful and that the differences between `port` and `marlin` come
from the fixes only.

Each implementation runs with two configurations: `default` (default settings) and `nosmear` (Poisson
smearing, electronic noise, threshold smearing and time smearing off). The remaining randomness is the
energy-loss fluctuations, so without smearing the cell IDs and times of the reconstructed hits must be
identical in all the implementations.

For every sub-detector, `analyse.py` reports the efficiency, the core width (IQR/1.349) and tails of the
residuals along the local u and v directions of the sensor, the cluster size, the cluster charge, the
time resolution, and a hit-by-hit comparison with Marlin (same cell ID, same time, same position, same
fired pixels). Detailed plots are made for one sub-detector and layer.

Both implementations use the MAIA_v0 geometry of the Marlin stack, because the endcap drivers changed in
the Key4hep stack. `ddsim` runs with `--enableDetailedShowerMode`: otherwise the momentum and path length
of the SimTrackerHits are not written to LCIO and both digitisers fall back on the MC particle momentum and
the straight crossing of the sensor.

## Running

Requirements: `git`, `python3`, and Docker or Apptainer with access to the
`ghcr.io/muoncollidersoft/mucoll-sim-ubuntu24` images.

```bash
test/validation/run_validation.sh
```

The steps can be run separately, e.g. to re-run the analysis only:

```bash
test/validation/run_validation.sh analyse
```

| Stage | Image | Content |
|---|---|---|
| `geometry` | both | Copy the MAIA_v0 compact files of the Marlin stack; export the sensor surfaces |
| `sim` | Marlin | Simulate the muon-gun sample; convert it to EDM4hep |
| `build` | both | Build the Marlin processor, the port and its variants |
| `digi` | both | Digitise with every implementation and configuration |
| `dump` | both | Dump the reconstructed hits, their SimTrackerHit and fired pixels to CSV |
| `analyse` | Key4hep | Write `summary.txt` and `validation_<subdet>_<layer>.png` |

The sample, the images and the implementations are set with environment variables, listed at the top of
`run_validation.sh`. For example:

```bash
WORKDIR=/tmp/validation NEVENTS=2000 THETA_MIN=30 THETA_MAX=150 \
ANALYSE_ARGS="--subdet ITBarrel --layer -1 --check" \
CONTAINER_RUNTIME=apptainer test/validation/run_validation.sh
```

Keep `WORKDIR` outside the source tree, or name it `validation_work` (ignored by git and by the copies of
the sources).

With `--check` (the default), the analysis fails if the port and Marlin differ in cell IDs or times
without smearing, or if `revert_all` is not compatible with Marlin (efficiency, or KS test with
p < 0.001 on the residuals, cluster size and charge).

## Results

1000 events of 10 muons, 1–100 GeV, θ between 30° and 150°, on layer 0 of the vertex barrel
(11624 SimTrackerHits), with default settings:

| | Efficiency | Residual core u / v (µm) | Beyond 15 µm u / v | Cluster size | Median charge |
|---|---|---|---|---|---|
| Marlin | 0.9820 | 4.80 / 4.79 | 0.2% / 0.3% | 3.05 | 17.3 keV |
| Gaudi port | 0.9822 | 5.17 / 5.13 | 0.5% / 0.7% | 2.97 | 18.9 keV |
| `CutOnDeltaRays` fix reverted | 0.9821 | 4.81 / 4.81 | 0.2% / 0.3% | 3.05 | 17.3 keV |
| All fixes reverted | 0.9821 | 4.78 / 4.79 | 0.2% / 0.4% | 3.05 | 17.3 keV |

- With the fixes reverted, the port is compatible with Marlin on every quantity (KS p ≥ 0.99). Without
  smearing, all the matched hits have the same cell ID and time; 87% have the same fired pixels, the
  rest differ because the energy-loss fluctuations use different random numbers.
- On the vertex barrel, the differences between the port and Marlin come from the `CutOnDeltaRays` fix:
  Marlin effectively ran with a delta-ray cut decreasing to about 1.7 keV instead of the configured
  30 keV, which suppressed the Landau tail of the energy loss. With the fix the charge is 9% higher, the
  clusters are 3% smaller and the residuals have larger tails.
- In the inner and outer tracker barrels, the ladder length fix moves the pixel grid by 10 µm in the
  local v direction, so that the pixel boundaries coincide with the edges of the 30.1 mm modules; the
  efficiency does not change. The fired pixels of the port and of Marlin therefore never have the same
  coordinates there (pixel overlap 0 unless `revert_all`).
- The same conclusions hold for all the tracker sub-detectors: with a 200-event sample, `revert_all` is
  compatible with Marlin everywhere, and the `CutOnDeltaRays` fix raises the median cluster charge by
  5–10% in every sub-detector.
- The Marlin processor crashes when the job terminates (`malloc_consolidate(): invalid chunk size`),
  after its output is complete.
