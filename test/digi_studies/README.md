# Studies of the tracker digitisation in MAIAConfig

These scripts compare the tracker hits that [MAIAConfig](https://github.com/MuonColliderSoft/MAIAConfig)
produces with its two digitisation flavours, and study the cluster shapes of the realistic one:

- **parametric**: `DDPlanarDigi` (k4Reco), configured in `TrackerDigi/tracking_{vertex,inner,outer}.py`,
  which smears the position and time of each `SimTrackerHit`;
- **realistic**: `MuonCVXDDigitiser` (this repository), configured in
  `TrackerDigi/realistic_{vertex,inner,outer}.py`, which transports the charge through the sensor and
  reconstructs each hit from a cluster of fired pixels.

The digitisers are built from MAIAConfig's own configuration functions, so the studies follow the
MAIAConfig checkout they are run with. Only the tracker digitisers run, not the full `digi_steer.py`. The
studies are not part of `ctest`.

## Running

Requirements: a MAIAConfig checkout, and Docker or Apptainer with access to the
`ghcr.io/muoncollidersoft/mucoll-sim-ubuntu24` images.

```bash
MAIACONFIG_SRC=/path/to/MAIAConfig test/digi_studies/run_digi_studies.sh
```

`MAIACONFIG_SRC` defaults to a `MAIAConfig` directory next to this repository. The stages can be run
separately, e.g. to redo the analysis only:

```bash
test/digi_studies/run_digi_studies.sh analyse
```

| Stage | Content |
|---|---|
| `build` | Build `MuonCVXDDigitiser`; `build/run` puts it in front of the one of the software stack |
| `geometry` | Export the sensor surfaces (`../validation/export_surfaces.py`) |
| `sim` | Simulate a muon-gun sample with `ddsim` |
| `digi` | Digitise the sample four times (see below) with `tracker_digi.py` |
| `dump` | Dump every `TrackerHitPlane` field, its `SimTrackerHit` and its fired pixels to CSV (`dump_hits.py`) |
| `analyse` | `compare.py` and `clusters.py`: summaries and plots |

The sample, the image and the container runtime are set with environment variables, listed at the top of
`run_digi_studies.sh`. Keep `WORKDIR` outside the source tree, or name it `digi_studies_work`.

The four digitisations:

| Name | Digitiser | Settings |
|---|---|---|
| `ddplanar` | `DDPlanarDigi` | MAIAConfig as is |
| `ddplanar_fixed` | `DDPlanarDigi` | with `CorrectTimeWindowForPropagation`, i.e. the time window cut applied to the time corrected for the time of flight, as intended by the configuration |
| `cvxd_long` | `MuonCVXDDigitiser` | MAIAConfig as is (long 25 µm × 312 µm tracker pixels), fired pixels stored |
| `cvxd_square` | `MuonCVXDDigitiser` | tracker `PixelSizeY` set to 25 µm, fired pixels stored |

## Outputs

`compare.py`: `DDPlanarDigi` vs `MuonCVXDDigitiser` as configured in MAIAConfig.

- `compare_summary.txt`: per sub-detector, the efficiency, the core width (IQR/1.349) and RMS of the
  residuals along the local u and v directions of the sensor, the stored uncertainties `du` and `dv`, the
  pulls (residual / `du`, `dv`, which is what ACTS uses), the time with respect to the `SimTrackerHit` with
  and without the time-of-flight correction, the energy, and which EDM fields each digitiser fills
  (`covMatrix`, `type`, `quality`, ...).
- `compare_resolution.png`: residuals and pulls along u and v.
- `compare_time_energy.png`: time convention, `EDep` / `EDep` of the `SimTrackerHit`, and efficiency
  against the time-of-flight corrected `SimTrackerHit` time.

`clusters.py`: cluster shapes of `MuonCVXDDigitiser` from the fired pixels, for square and long tracker
pixels, and resolutions against `DDPlanarDigi` (with the time window corrected).

- `clusters_summary.txt`: per sub-detector, the efficiency, cluster size and extent along u and v,
  resolutions, uncertainties, pulls, and the spread of the pixel times in a cluster.
- `clusters_shapes.png`: cluster size and extent, and the mean extent against the incidence angle with the
  geometric expectation: thickness × |tan α| / pitch + 1 along v, thickness × |tan α − tan θ_L| / pitch + 1
  along u, where the Lorentz drift (`TanLorentz` = 0.8) shears the charge.
- `clusters_resolution_v.png`: residuals and pulls along v, and the charge per fired pixel.
- `clusters_examples.png`: example clusters in bins of incidence angle, with the charge of every pixel,
  drawn with the physical aspect ratio of the pixels.

## Fired pixels in the EDM

With `StoreFiredPixels` (`--storeFiredPixels` in MAIAConfig), `MuonCVXDDigitiser` writes for every hit its
fired pixels, as `SimTrackerHit`s in `<prefix>Pixels`, and links from the hit to its pixels in
`<prefix>RawHitRelations` (weight 1/N). A pixel has its position in pixel units in the local (u, v) frame of
the sensor (z = 0), its digitised charge in electrons as `EDep`, its own smeared and digitised time, and the
MC particle and flags of the `SimTrackerHit`.

## Results

300 events of 10 muons, 1–100 GeV, all θ, MAIA_v0, `post_3_1_test` software stack:

- `DDPlanarDigi` as MAIAConfig configures it keeps only about a third of the vertex barrel hits and almost
  none of the others: since k4Reco #65 the time window cut is corrected for the time of flight only with
  `CorrectTimeWindowForPropagation`, which MAIAConfig does not set, so the window is applied to the raw
  hit times. With it, 89–98% of the `SimTrackerHit`s are kept, against 94–98% for `MuonCVXDDigitiser`.
- Both digitisers write the same collections with the same cell IDs and U/V. `DDPlanarDigi` fills
  `covMatrix`, stores the time corrected for the time of flight (t − r/c, smeared) and the `SimTrackerHit`
  energy. `MuonCVXDDigitiser` leaves `covMatrix` empty, stores the absolute time (smeared and digitised in
  9.8 ps steps up to 10 ns), applies no time window, and stores the digitised cluster charge.
- Vertex detector: both give about 5 µm along u and v. `MuonCVXDDigitiser` stores du = dv = pitch / √12 =
  7.2 µm, so its pulls have a width of 0.6–0.8.
- Trackers, with long pixels: along v, 81–96 µm against 88–91 µm for `DDPlanarDigi`, pulls 0.9–1.07; along
  u, 6.3–6.5 µm against about 7 µm, pulls about 0.9. With square 25 µm pixels the trackers were about 20
  times better along v than the parametric digitisation assumes.
- Cluster extents follow the geometric expectation along u and v. The pixel times in a cluster spread by
  100–300 ps, as each pixel time is smeared independently.
