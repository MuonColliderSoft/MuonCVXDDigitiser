# Tests

Regression tests of `MuonCVXDDigitiser` in Marlin, run with `ctest` in an environment with Marlin,
LCIO and DD4hep and `K4GEO` set: the Key4hep stack (as in CI), or e.g.
`ghcr.io/muoncollidersoft/mucoll-sim-ubuntu24:v2.11`:

```bash
cmake -S . -B build && cmake --build build -j8
ctest --test-dir build -j4 --output-on-failure
```

## Samples

Two muon-gun samples in the MAIA_v0 geometry (10 μ⁻ per event, 1–100 GeV, θ 10°–170°) are attached
to the [`lcio-test-inputs`](https://github.com/MuonColliderSoft/MuonCVXDDigitiser/releases/tag/lcio-test-inputs)
release and downloaded by the first test, `fetch_test_data`, which checks their SHA-256:

| File | Events | Content |
|---|---|---|
| `mu_gun_detailed.slcio` | 10 | `ddsim --enableDetailedShowerMode`: hits with momentum and path length |
| `mu_gun_nodetail.slcio` | 5 | Default `ddsim`: hits without momentum and path length, to test the fallbacks |

To work offline, configure with `-DMUONCVXD_TEST_DATA_DIR=<folder with the samples>`.

## Tests

| Test | What it does |
|---|---|
| `digi_<name>` | Runs the digitiser in Marlin with a configuration of `digi_common.py` and dumps the hits to CSV |
| `check_<name>` | Sanity requirements (finite values, fired pixels, efficiency, hits near their SimTrackerHit, residuals), then the summary statistics against `reference/<name>.json` |
| `compare_reproducible` | The same configuration run twice gives identical hits |

Only `libMarlinDD4hep` and this package are loaded: other plugins of the software stack crash when
the job terminates, which would hide failures.

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
