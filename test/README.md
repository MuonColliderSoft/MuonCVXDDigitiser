# Tests

Unit tests and regression tests of the `MuonCVXDDigitiser` and `MuonCVXDRealDigitiser` Gaudi
algorithms, run with `ctest` in a Key4hep environment with `K4GEO` set: the Key4hep stack (as in
CI), or e.g. `ghcr.io/muoncollidersoft/mucoll-sim-ubuntu24:main` (start it with
`--entrypoint bash`: its default entrypoint ignores the command):

```bash
cmake -S . -B build && cmake --build build -j8
ctest --test-dir build -j4 --output-on-failure
```

## Samples

Muon-gun samples in the MAIA_v0 geometry (10 μ⁻ per event, 1–100 GeV, θ 10°–170°) are attached
to the [`edm4hep-test-inputs`](https://github.com/MuonColliderSoft/MuonCVXDDigitiser/releases/tag/edm4hep-test-inputs)
release and downloaded by the first test, `fetch_test_data`, which checks their SHA-256:

| File | Events | Content |
|---|---|---|
| `mu_gun_detailed.edm4hep.root` | 10 | `ddsim --enableDetailedShowerMode`: hits with momentum and path length |
| `mu_gun_nodetail.edm4hep.root` | 5 | Default `ddsim`: hits without momentum and path length, to test the fallbacks |
| `mu_gun_scaled.edm4hep.root` | 10 | `mu_gun_detailed` with the Geant4 deposit doubled and the path length halved |

They are the LCIO samples of the [`lcio-test-inputs`](https://github.com/MuonColliderSoft/MuonCVXDDigitiser/releases/tag/lcio-test-inputs)
release, on which the references of the Marlin processors were made, converted with
`lcio2edm4hep`.

To work offline, configure with `-DMUONCVXD_TEST_DATA_DIR=<folder with the samples>`.

## Tests

| Test | What it does |
|---|---|
| `test_charge_bins`, `test_sensors`, `test_shared_infrastructure` | Unit tests: charge digitization bins, sensor models, geometry loading, cell IDs and energy-loss fluctuations |
| `digi_<name>` | Writes a `k4run` options file with a configuration of `digi_common.py` (`make_options.py`), runs it and dumps the hits to CSV (`dump_hits.py`) |
| `check_<name>` | Sanity requirements (finite values, fired pixels, efficiency, hits near their SimTrackerHit, residuals), then the summary statistics against `reference/<name>.json` |
| `compare_reproducible` | The same configuration run twice gives identical hits |
| `compare_default_mt_edm4hep`, `compare_real_mt_edm4hep` | Run with 4 threads, every hit, fired pixel and link is identical to the single-threaded run |
| `check_default_edm4hep`, `check_nopixels_edm4hep`, `check_real_<barrel>_edm4hep` | Consistency of the EDM4hep output: links to SimTrackerHits on the same sensor, raw hit links to every fired pixel with weights summing to 1, no pixels unless `StoreFiredPixels` |
| `compare_timing_overrides` | The realistic time smearing with every term overridden gives the same hits as the constant model with the same sigma |
| `compare_resim_ignores_geant4` | On a copy of the sample with the Geant4 deposit doubled and the path length halved (`make_scaled_sample`), the re-simulated ionisation trail gives identical hits |
| `compare_geant4_charge_follows_deposit` | On the same copy, the Geant4 ionisation trail gives twice the charge |
| `check_time_bins` | Two digitisers in one job with different time digitizations: every fired-pixel time is on the bins of its own digitiser |
| `digi_real`, `digi_real_rd53a` | `MuonCVXDRealDigitiser` on the three barrels (it does not support endcaps), with the trivial and the RD53A sensor models. Its hits are linked to all their SimTrackerHits, so `mean_pixels` is the number of SimTrackerHits per cluster |
| `compare_real_reproducible` | The same for `MuonCVXDRealDigitiser` |
| `options_<algorithm>` | The example options files in `MuonCVXDDigitiser/options` run on the sample |

## Two levels of regression checks

- **References** (`reference/*.json`, `check_<name>`): summary statistics per sub-detector with
  tolerances of a few percent, so they hold on any platform. Coarse: small changes can stay within
  the tolerances.
CI (`.github/workflows/key4hep-build.yaml`) builds with `-Werror` and runs the tests against the
Key4hep release (required) and nightlies (reported, not blocking), from CVMFS.

- **Comparison with the base branch** (CI, `compare_to_base.sh`): on a pull request, CI also builds
  the base branch and compares the outputs hit by hit. The tables are in the job summary, and the
  result is an annotation on the run:

  | Outputs | Without `changes-output` | With `changes-output` |
  |---|---|---|
  | identical | passes | passes, with a warning to remove the label |
  | different | fails (error) | passes, with a notice |
  | comparison failed (build or digitisation) | fails (error) | fails (error) |

  The label applies to the whole pull request: when only some of its commits are meant to change
  the output, check the others locally with `compare_to_base.sh`.

A pull request that changes the output on purpose regenerates the references, so the change is
visible in its diff, and gets the `changes-output` label:

```bash
ctest --test-dir build -j4
cmake --build build --target update_test_references
```
