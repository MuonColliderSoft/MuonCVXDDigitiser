# MuonCVXDDigitiser

Realistic digitisation of pixelated silicon sensors for the Muon Collider, as
[Gaudi](https://gitlab.cern.ch/gaudi/Gaudi) algorithms for the
[Key4hep](https://key4hep.github.io/key4hep-doc/) framework
([k4FWCore](https://github.com/key4hep/k4FWCore), [EDM4hep](https://github.com/key4hep/EDM4hep)).

Starting from Geant4 `SimTrackerHit`s, the charge released along the track in the sensor is
drifted to the readout plane (Lorentz angle, diffusion) and shared among pixels. The pixel
response is then simulated and a `TrackerHitPlane` is reconstructed from each cluster of
fired pixels. The procedure follows the
[CMS pixel digitisation](https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePixelDigitization).

The package provides two algorithms:

| Algorithm | Sub-detectors | Model |
|---|---|---|
| `MuonCVXDDigitiser` | Vertex, inner and outer tracker; barrels and endcaps | Each `SimTrackerHit` is digitised on its own and gives at most one `TrackerHitPlane`. Poisson and electronic noise, smeared threshold, charge and time discretisation, and a time resolution model. |
| `MuonCVXDRealDigitiser` | Barrel only; practical for the vertex barrel | All the `SimTrackerHit`s of a ladder are fed in time order to a model of the readout chip through a sliding time window. Pixels are clustered on each sensor, so hits from different particles can merge. |

This is the port to Key4hep of the Marlin processors of
[spg-berkeleylab/MuonCVXDDigitiser](https://github.com/spg-berkeleylab/MuonCVXDDigitiser);
see [Changes with respect to the Marlin processors](#changes-with-respect-to-the-marlin-processors).

## Building

The package needs Gaudi, k4FWCore, EDM4hep, podio, DD4hep (with DDRec), ROOT, CLHEP and GSL,
all provided by the Key4hep and Muon Collider software stacks. In an environment where the
stack is set up (e.g. the `ghcr.io/muoncollidersoft/mucoll-sim-ubuntu24` images, after
`source /opt/setup_mucoll.sh`):

```bash
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=install
cmake --build build -j
cmake --install build
```

Then put the installation in front of the stack in the environment:

```bash
export LD_LIBRARY_PATH=$PWD/install/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$PWD/install/python:$PYTHONPATH
export GAUDI_PLUGIN_PATH=$PWD/install/lib:$GAUDI_PLUGIN_PATH
```

Without installing, `build/run <command>` runs a command in the environment of the build tree.

The API documentation is built with Doxygen by configuring with `-DINSTALL_DOC=ON`.

## Testing

```bash
cd build
ctest -j4 --output-on-failure
```

The tests use the MAIA geometry from `$K4GEO` and take a few minutes: a small muon-gun sample is
simulated with `ddsim`, then digitised with both algorithms. The outputs are checked (efficiency,
residuals, links, cell IDs) and runs with 1 and 4 threads must give identical results.
Unit tests cover the geometry loading, the cell ID coding, the energy-loss fluctuation models and
the pixel clustering.

## Running

Example options files for the MAIA detector are in `MuonCVXDDigitiser/options`:

```bash
# All vertex and tracker sub-detectors with MuonCVXDDigitiser
k4run MuonCVXDDigitiser/options/runMuonCVXDDigitiser.py \
      --IOSvc.Input sim.edm4hep.root --IOSvc.Output digi.edm4hep.root [--threads N]

# The vertex barrel with MuonCVXDRealDigitiser
k4run MuonCVXDDigitiser/options/runMuonCVXDRealDigitiser.py \
      --IOSvc.Input sim.edm4hep.root --IOSvc.Output digi_real.edm4hep.root \
      [--sensor-type {0,1}] [--stats-file stats.root] [--threads N]
```

Both algorithms need the `GeoSvc`, to read the geometry, and the `UniqueIDGenSvc`, to seed the
random numbers of each event from the event and run numbers. The output is therefore
reproducible and does not depend on the number of threads.

A minimal configuration:

```python
from Configurables import MuonCVXDDigitiser

digitiser = MuonCVXDDigitiser(
    "VXDBarrelDigitiser",
    SubDetectorName="VertexBarrel",
    LayerIDs=[0, 1, 2, 4, 6],
    CollectionName=["VertexBarrelCollection"],
    EventHeader=["EventHeader"],
    OutputCollectionName=["VXDBarrelHits"],
    RelationColName=["VXDBarrelHitsRelations"],
    SimHitLocCollectionName=["VXDBarrelPixels"],
    RawHitsLinkColName=["VXDBarrelRawHitRelations"],
)
```

`LayerIDs` lists the values of the `layer` field of the cell IDs, in the order of the layers of the
geometry. For MAIA_v0 they are `[0, 1, 2, 4, 6]` for the vertex barrel, `[0, ..., 7]` for the vertex
endcap, `[0, 1, 2]` for the inner and outer tracker barrels, `[0, ..., 6]` for the inner tracker
endcap and `[0, 1, 2, 3]` for the outer tracker endcap.

## MuonCVXDDigitiser

### Collections

| Property | Direction | Type | Content |
|---|---|---|---|
| `CollectionName` | input | `SimTrackerHitCollection` | Simulated hits of the sub-detector |
| `EventHeader` | input | `EventHeaderCollection` | Event and run numbers, to seed the random numbers |
| `OutputCollectionName` | output | `TrackerHitPlaneCollection` | Reconstructed hits |
| `RelationColName` | output | `TrackerHitSimTrackerHitLinkCollection` | Link from each reconstructed hit to its `SimTrackerHit`, weight 1 |
| `SimHitLocCollectionName` | output | `SimTrackerHitCollection` | Fired pixels, filled only with `StoreFiredPixels`: position in pixel units in the local frame of the sensor (x, y; z = 0), `EDep` in electrons, and the time, MC particle and flags of the original `SimTrackerHit` |
| `RawHitsLinkColName` | output | `TrackerHitSimTrackerHitLinkCollection` | Link from each reconstructed hit to its fired pixels, weight 1/(number of pixels), filled only with `StoreFiredPixels` |

The reconstructed hits have the cell ID of their `SimTrackerHit`, `EDep` in GeV, the time in ns,
`du` = `PixelSizeX`/√12 and `dv` = `PixelSizeY`/√12. `u` and `v` hold the (θ, φ) of the directions
of the sensor surface.

### Properties

| Property | Type | Default | Description |
|---|---|---|---|
| `SubDetectorName` | string | `"VertexBarrel"` | Name of the sub-detector in the geometry. Must contain `Barrel` or `Endcap`, and `Vertex`, `InnerTracker` or `OuterTracker` |
| `LayerIDs` | vector<int> | `[]` | Cell ID layer IDs of the geometry layers (required) |
| `EncodingStringParameterName` | string | `"GlobalTrackerReadoutID"` | DD4hep constant with the cell ID encoding |
| `GeoSvcName` | string | `"GeoSvc"` | Name of the GeoSvc instance |
| `ZSegmented` | int | `-1` | Sensor segmentation along z, barrel layers only: -1 = auto (on for the vertex barrel), 0 = off, 1 = on |
| `PixelSizeX` | double | `0.025` | Pixel size across the ladder (mm) |
| `PixelSizeY` | double | `0.025` | Pixel size along the ladder (mm) |
| `TanLorentz` | double | `0.8` | Tangent of the Lorentz angle |
| `TanLorentzY` | double | `0.` | Tangent of the Lorentz angle along y |
| `DiffusionCoefficient` | double | `0.07` | Diffusion coefficient, √(2D/μ/V) |
| `ElectronsPerKeV` | double | `270.3` | Electrons per keV of deposited energy |
| `CutOnDeltaRays` | double | `0.030` | Cut on delta-ray energy (MeV) |
| `SegmentLength` | double | `0.005` | Length of the segments of the ionisation trail (mm) |
| `MaxTrackLength` | double | `10.0` | Maximum track length in the sensor (mm) |
| `MaxEnergyDelta` | double | `100.0` | Maximum difference between the Geant4 deposit and the sampled charge (electrons) |
| `ResimulateIonisation` | bool | `false` | `false`: path length and deposited energy from Geant4; `true`: re-simulated from the local direction and `EnergyLoss` |
| `EnergyLoss` | double | `280.0` | Mean energy loss (keV/mm), used with `ResimulateIonisation` or for hits without deposit |
| `DoMultipleScattering` | bool | `false` | Multiple scattering of the track in the sensor |
| `MSSliceThickness` | double | `0.005` | Thickness of the slices for the multiple scattering (mm) |
| `PoissonSmearing` | bool | `true` | Poisson smearing of the pixel charge |
| `ElectronicEffects` | bool | `true` | Gaussian electronic noise on the pixel charge |
| `ElectronicNoise` | double | `80.` | Electronic noise (electrons) |
| `Threshold` | double | `500.` | Pixel threshold (electrons) |
| `ThresholdSmearSigma` | int | `25` | Sigma of the Gaussian smearing of the threshold of each pixel (electrons) |
| `DigitizeCharge` | bool | `true` | Discretisation of the pixel charge |
| `ChargeDigitizeNumBits` | int | `4` | Number of bits of the charge discretisation |
| `ChargeDigitizeBinning` | int | `1` | Charge binning: 0 = uniform between `Threshold` and `ChargeMaximum`, 1 = variable, from tables for 3, 4, 5, 6 and 8 bits (dedicated 4-bit tables for vertex sensors of 50, 75, 100, 200 and 400 µm; bins twice as wide for the trackers) |
| `ChargeMaximum` | double | `15000.` | Dynamic range of the pixel charge (electrons) |
| `TimeSmearingModel` | int | `2` | Time resolution: 0 = none, 1 = constant `TimeSmearingSigma`, 2 = realistic, from the sensor thickness |
| `TimeSmearingSigma` | double | `0.05` | Time resolution for `TimeSmearingModel` = 1 (ns) |
| `TRise`, `SigmaLandau`, `SigmaTimewalk`, `SigmaJitter`, `SigmaTDC`, `SigmaClock` | double | `-1.0` | Overrides of the terms of the realistic time resolution (ns); negative = derived from the sensor thickness |
| `DigitizeTime` | bool | `true` | Discretisation of the pixel time |
| `TimeDigitizeNumBits` | int | `10` | Number of bits of the time discretisation |
| `TimeDigitizeBinning` | int | `0` | Time binning: 0 = uniform up to `TimeMaximum` (the only scheme) |
| `TimeMaximum` | double | `10.0` | Dynamic range of the pixel time (ns) |
| `StoreFiredPixels` | bool | `false` | Fill the fired pixel and raw hit link collections |

## MuonCVXDRealDigitiser

### Collections

| Property | Direction | Type | Content |
|---|---|---|---|
| `CollectionName` | input | `SimTrackerHitCollection` | Simulated hits of the sub-detector |
| `EventHeader` | input | `EventHeaderCollection` | Event and run numbers, to seed the random numbers |
| `OutputCollectionName` | output | `TrackerHitPlaneCollection` | Reconstructed hits, one per cluster |
| `RelationColName` | output | `TrackerHitSimTrackerHitLinkCollection` | Links from each reconstructed hit to all the `SimTrackerHit`s contributing to its cluster, weight 1 |

### Properties

| Property | Type | Default | Description |
|---|---|---|---|
| `SubDetectorName` | string | `"VertexBarrel"` | Name of the barrel sub-detector in the geometry |
| `LayerIDs` | vector<int> | `[]` | Cell ID layer IDs of the geometry layers; empty = 0, 1, 2, ... |
| `EncodingStringParameterName` | string | `"GlobalTrackerReadoutID"` | DD4hep constant with the cell ID encoding |
| `GeoSvcName` | string | `"GeoSvc"` | Name of the GeoSvc instance |
| `ZSegmented` | int | `-1` | Sensor segmentation along z: -1 = auto (on for the vertex barrel), 0 = off, 1 = on |
| `SensorType` | int | `1` | Sensor model: 0 = RD53A chip (charge from time over threshold), 1 = trivial (charge collected in the time window) |
| `WindowSize` | float | `25.` | Clock period, i.e. width of the time window (ns) |
| `RD53Aslope` | float | `0.1` | Discharge slope of the RD53A front end (electrons/ns) |
| `PixelSizeX` | double | `0.025` | Pixel size across the ladder (mm) |
| `PixelSizeY` | double | `0.025` | Pixel size along the ladder (mm) |
| `Threshold` | double | `200.` | Pixel threshold (electrons) |
| `TanLorentz` | double | `0.8` | Tangent of the Lorentz angle |
| `TanLorentzY` | double | `0.` | Tangent of the Lorentz angle along y |
| `DiffusionCoefficient` | double | `0.07` | Diffusion coefficient, √(2D/μ/V) |
| `ElectronsPerKeV` | double | `270.3` | Electrons per keV of deposited energy |
| `CutOnDeltaRays` | double | `0.030` | Cut on delta-ray energy (MeV) |
| `SegmentLength` | double | `0.005` | Length of the segments of the ionisation trail (mm) |
| `EnergyLoss` | double | `280.0` | Mean energy loss (keV/mm) |
| `MaxEnergyDelta` | double | `100.0` | Maximum difference between the Geant4 deposit and the sampled charge (electrons) |
| `MaxTrackLength` | double | `10.0` | Maximum track length in the sensor (mm) |
| `CreateStats` | bool | `false` | Fill cluster size, cluster extent, cluster energy and hit offset histograms, separately for signal and beam-induced background (clusters with only overlay hits). Write them with a histogram sink, e.g. `Gaudi::Histograming::Sink::Root` |

The pixel matrix of a whole ladder is kept in memory: about 11 MB for a MAIA vertex barrel
ladder, but more than 100 MB for a tracker barrel ladder.

## Changes with respect to the Marlin processors

The physics follows spg-berkeleylab/MuonCVXDDigitiser `master` (31bb9e2). The interface of
`MuonCVXDDigitiser` is the one already used by
[MAIAConfig](https://github.com/MuonColliderSoft/MAIAConfig).

### Interface

- LCIO `LCRelation`s are replaced by EDM4hep link collections. `TrackerHitPlane` has no raw hits in
  EDM4hep: the fired pixels of `MuonCVXDDigitiser` are linked through `RawHitsLinkColName`.
- The fired pixel collection is named by `SimHitLocCollectionName` instead of the fixed
  `VBPixels`, `VEPixels`, `IBPixels`, `IEPixels`, `OBPixels` and `OEPixels`.
- Output cell IDs use the encoding of the geometry, with all 64 bits, instead of the LCIO
  `LCTrackerCellID` encoding.
- On/off parameters are `bool` properties.
- `MuonCVXDRealDigitiser`: `LayerIDs` is new; `StatisticsFilename` is replaced by `CreateStats` and
  the histogram sink; `PoissonSmearing`, `ElectronicEffects`, `ElectronicNoise` and
  `StoreFiredPixels` are removed, as they had no effect.
- The compile-time `ZSEGMENTED` option and OpenMP are gone: segmentation is set by `ZSegmented`,
  and parallelism comes from the Gaudi multithreaded event processing.

### Fixes that change the results

`MuonCVXDDigitiser`:
- The ladder length of tracker barrels, which have no sensor length in the geometry, is converted
  from cm to mm: it was 10 times too short.
- Each pixel has its own threshold, smeared around `Threshold`, instead of a smeared threshold
  accumulated from pixel to pixel along the cluster.
- `CutOnDeltaRays` stays constant: the fluctuation model modified it in place for later segments,
  hits and events.
- Hits whose layer ID is not in `LayerIDs` are skipped before any geometry lookup.
- The time discretisation of every instance uses its own `TimeDigitizeNumBits` and `TimeMaximum`
  (they were shared by all instances in the job).
- The z of the fired pixels is 0 instead of uninitialised.

`MuonCVXDRealDigitiser`:
- Layers are identified through `LayerIDs`: the geometry layer index was used as the cell ID layer,
  which mixes up the geometry of the MAIA vertex barrel layers 4 and 6.
- The cell IDs of the clusters are encoded with the `system` field of the geometry encoding.
- Pixel clustering no longer reads outside the pixel matrix nor merges pixels across rows or across
  the sensors of a ladder.
- Layers without z segmentation (e.g. tracker barrels) use one sensor per ladder instead of dividing
  by zero.
- Hits and clusters without a surface are skipped.
- Statistics: the cluster extents are computed correctly for single-hit clusters, the RD53A sensor
  model fills the cluster size, and the BIB energy histogram has its own name.

## Package layout

```
MuonCVXDDigitiser/
  components/   Gaudi algorithms
  include/ src/ Framework-independent code: geometry, cell IDs, random numbers,
                energy-loss fluctuations, sensor models (library MuonCVXDDigitiserCore)
  options/      Example k4run options
test/           Unit tests and ctest checks
cmake/          CMake helpers and package configuration
doc/            Doxygen configuration
```
