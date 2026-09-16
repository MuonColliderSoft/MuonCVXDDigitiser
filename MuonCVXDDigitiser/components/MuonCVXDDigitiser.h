#ifndef MuonCVXDDigitiser_h
#define MuonCVXDDigitiser_h 1

#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include "Gaudi/Accumulators.h"
#include "Gaudi/Property.h"
#include "k4FWCore/Transformer.h"
#include "k4Interface/IGeoSvc.h"
#include "k4Interface/IUniqueIDGenSvc.h"

#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/TrackerHitPlaneCollection.h"
#include "edm4hep/TrackerHitSimTrackerHitLinkCollection.h"

#include "DigitiserTypes.h"
#include "LayerGeometry.h"
#include "MyG4UniversalFluctuationForSi.h"
#include "TrackerCellID.h"

class EventRandom;

/** Digitizer for Simulated Hits in the Vertex and Tracker Detectors. <br>
 * Digitization follows the procedure adopted in the CMS software package.
 * See https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePixelDigitization
 *
 * For every SimTrackerHit the ionisation trail inside the sensor is split into segments,
 * the charge is drifted (Lorentz angle, diffusion) to the readout plane and shared among
 * pixels. Pixel charges are smeared (Poisson, electronic noise), thresholded and
 * discretised, pixel times are smeared and discretised, and a TrackerHitPlane is
 * reconstructed from the cluster.
 *
 * Inputs:
 * - CollectionName: SimTrackerHits of the sub-detector
 * - EventHeader: used to seed the random numbers for each event
 *
 * Outputs:
 * - OutputCollectionName: reconstructed TrackerHitPlanes
 * - RelationColName: links from each reconstructed hit to its SimTrackerHit (weight 1)
 * - SimHitLocCollectionName: fired pixels, as SimTrackerHits with the position in pixel units in
 *   the local sensor frame and EDep in electrons. Only filled if StoreFiredPixels is set.
 * - RawHitsLinkColName: links from each reconstructed hit to its fired pixels (weight 1/N pixels).
 *   Only filled if StoreFiredPixels is set.
 *
 * All the parameters are documented in the property declarations below.
 */
struct MuonCVXDDigitiser final
    : k4FWCore::MultiTransformer<std::tuple<edm4hep::SimTrackerHitCollection,
                                            edm4hep::TrackerHitPlaneCollection,
                                            edm4hep::TrackerHitSimTrackerHitLinkCollection,
                                            edm4hep::TrackerHitSimTrackerHitLinkCollection>(
          const edm4hep::SimTrackerHitCollection&, const edm4hep::EventHeaderCollection&)>
{
    MuonCVXDDigitiser(const std::string& name, ISvcLocator* svcLoc);

    StatusCode initialize() override;

    std::tuple<edm4hep::SimTrackerHitCollection,
               edm4hep::TrackerHitPlaneCollection,
               edm4hep::TrackerHitSimTrackerHitLinkCollection,
               edm4hep::TrackerHitSimTrackerHitLinkCollection>
    operator()(const edm4hep::SimTrackerHitCollection& simTrackerHits,
               const edm4hep::EventHeaderCollection& headers) const override;

private:
    /// A pixel fired by the current SimTrackerHit. Position in local sensor coordinates (mm).
    struct PixelHit
    {
        double x;
        double y;
        double z;
        float charge; // electrons
        float time;   // ns
    };
    /// Pixels keyed on the pixel index: the ordering fixes the order in which random numbers are drawn
    typedef std::map<int, PixelHit> PixelHitMap;

    /// Reconstructed cluster, in local sensor coordinates
    struct ClusterHit
    {
        double position[3];
        float eDep; // DD4hep energy units
        float time;
    };

    /// Digitisation state of the SimTrackerHit being processed
    struct HitState
    {
        int currentLayer{0};
        int currentLadder{0};
        int numberOfSegments{0};
        double currentParticleMass{0};
        double currentParticleMomentum{0};
        double currentPhi{0};
        double eSum{0};
        double segmentDepth{0};
        // Pixel matrix of the sensor: [-sensorHalfLengthU, sensorHalfLengthU] x [-sensorHalfLengthV, sensorHalfLengthV]
        double sensorHalfLengthU{0};
        double sensorHalfLengthV{0};
        double currentLocalPosition[3]{};
        double currentEntryPoint[3]{};
        double currentExitPoint[3]{};
        IonisationPointVec ionisationPoints;
        SignalPointVec signalPoints;
    };

    Gaudi::Property<std::string> m_subDetName{this, "SubDetectorName", "VertexBarrel", "Name of Vertex detector"};
    Gaudi::Property<std::vector<int>> m_layerIDs{this, "LayerIDs", {}, "ID of layers of subdetector"};
    Gaudi::Property<std::string> m_encodingStringVariable{this, "EncodingStringParameterName", "GlobalTrackerReadoutID",
        "The name of the DD4hep constant that contains the cell ID encoding string for the sub-detector"};
    Gaudi::Property<std::string> m_geoSvcName{this, "GeoSvcName", "GeoSvc", "The name of the GeoSvc instance"};

    Gaudi::Property<double> m_tanLorentzAngleX{this, "TanLorentz", 0.8, "Tangent of Lorentz Angle"};
    Gaudi::Property<double> m_tanLorentzAngleY{this, "TanLorentzY", 0., "Tangent of Lorentz Angle along Y"};
    Gaudi::Property<double> m_cutOnDeltaRays{this, "CutOnDeltaRays", 0.030, "Cut on delta-ray energy (MeV)"};
    // For diffusion-coeffieint calculation, see e.g. https://www.slac.stanford.edu/econf/C060717/papers/L008.PDF
    // or directly Eq. 13 of https://cds.cern.ch/record/2161627/files/ieee-tns-07272141.pdf
    // diffusionCoefficient = sqrt(2*D / mu / V), where
    //  - D = 12 cm^2/s // diffusion constant
    //  - mu = 450 cm^2/s/V // mobility
    //  - V = 10-30 V // expected depletion voltage
    //  => _diffusionCoefficient = 0.04-0.07
    Gaudi::Property<double> m_diffusionCoefficient{this, "DiffusionCoefficient", 0.07, "Diffusion coefficient, sqrt(D / mu / V)."};
    Gaudi::Property<double> m_pixelSizeX{this, "PixelSizeX", 0.025, "Pixel size along the local u direction of the sensor (mm)"};
    Gaudi::Property<double> m_pixelSizeY{this, "PixelSizeY", 0.025, "Pixel size along the local v direction of the sensor (mm)"};
    Gaudi::Property<double> m_electronsPerKeV{this, "ElectronsPerKeV", 270.3, "Electrons per keV"};
    Gaudi::Property<double> m_threshold{this, "Threshold", 500., "Cell Threshold in electrons"};
    Gaudi::Property<double> m_chargeMax{this, "ChargeMaximum", 15000., "Cell dynamic range in electrons"};
    Gaudi::Property<double> m_segmentLength{this, "SegmentLength", 0.005, "Segment Length in mm"};
    Gaudi::Property<bool> m_PoissonSmearing{this, "PoissonSmearing", true, "Apply Poisson smearing of electrons collected on pixels"};
    Gaudi::Property<int> m_thresholdSmearSigma{this, "ThresholdSmearSigma", 25, "sigma of Gaussian used in threshold smearing, in electrons"};
    Gaudi::Property<bool> m_DigitizeCharge{this, "DigitizeCharge", true, "Flag to enable Digitization of the charge collected on pixels"};
    Gaudi::Property<int> m_ChargeDigitizeNumBits{this, "ChargeDigitizeNumBits", 4, "Number of bits used to determine bins for charge discretization"};
    Gaudi::Property<int> m_ChargeDigitizeBinning{this, "ChargeDigitizeBinning", 1, "Binning scheme used for charge discretization"};
    Gaudi::Property<bool> m_DigitizeTime{this, "DigitizeTime", true, "Flag to enable digitization of timing information."};
    Gaudi::Property<int> m_TimeDigitizeNumBits{this, "TimeDigitizeNumBits", 10, "Number of bits used to determine bins for time discretization"};
    Gaudi::Property<int> m_TimeDigitizeBinning{this, "TimeDigitizeBinning", 0, "Binning scheme used for time discretization"};
    Gaudi::Property<double> m_timeMax{this, "TimeMaximum", 10.000, "Cell dynamic range for timing measurement [ns]"};
    Gaudi::Property<double> m_timeSmearingSigma{this, "TimeSmearingSigma", 0.05,
        "Constant intrinsic time measurement resolution (ns), used when TimeSmearingModel = 1."};
    Gaudi::Property<int> m_timeSmearingModel{this, "TimeSmearingModel", 2,
        "Time smearing model: 0 = none, 1 = constant sigma (TimeSmearingSigma), 2 = realistic, derived from sensor thickness."};
    Gaudi::Property<bool> m_electronicEffects{this, "ElectronicEffects", true, "Apply Electronic Effects"};
    Gaudi::Property<double> m_electronicNoise{this, "ElectronicNoise", 80., "electronic noise in electrons"};
    Gaudi::Property<bool> m_produceFullPattern{this, "StoreFiredPixels", false, "Store fired pixels"};
    Gaudi::Property<double> m_energyLoss{this, "EnergyLoss", 280.0, "Energy Loss keV/mm"};
    Gaudi::Property<double> m_deltaEne{this, "MaxEnergyDelta", 100.0,
        "Max delta in energy between G4 prediction and random sampling for each hit in electrons"};
    Gaudi::Property<double> m_maxTrkLen{this, "MaxTrackLength", 10.0, "Maximum values for track length (in mm)"};
    Gaudi::Property<bool> m_resimulateIonisation{this, "ResimulateIonisation", false,
        "Source of the ionisation trail: false = take the path length and the deposited energy from the Geant4 hit, "
        "true = re-simulate them from the local direction and the EnergyLoss parametrisation."};
    Gaudi::Property<bool> m_doMultipleScattering{this, "DoMultipleScattering", false,
        "Flag to enable multiple scattering of the track inside the sensor."};
    Gaudi::Property<double> m_msSliceThickness{this, "MSSliceThickness", 0.005,
        "Slice thickness (mm) used to step the track through the sensor when applying multiple scattering."};
    Gaudi::Property<double> m_t_riseOverride{this, "TRise", -1.0, "Optional override for t_rise (ns). Negative means use default."};
    Gaudi::Property<double> m_sigma_landauOverride{this, "SigmaLandau", -1.0, "Optional override for sigma_landau (ns). Negative means use default."};
    Gaudi::Property<double> m_sigma_timewalkOverride{this, "SigmaTimewalk", -1.0, "Optional override for sigma_timewalk (ns). Negative means use default."};
    Gaudi::Property<double> m_sigma_jitterOverride{this, "SigmaJitter", -1.0, "Optional override for sigma_jitter (ns). Negative means use default."};
    Gaudi::Property<double> m_sigma_TDCOverride{this, "SigmaTDC", -1.0, "Optional override for sigma_TDC (ns). Negative means use default."};
    Gaudi::Property<double> m_sigma_clockOverride{this, "SigmaClock", -1.0, "Optional override for sigma_clock (ns). Negative means use default."};

    SmartIF<IGeoSvc> m_geoSvc;
    SmartIF<IUniqueIDGenSvc> m_uidSvc;

    // Set in initialize() and read-only afterwards
    LayerGeometry m_geo;
    std::unique_ptr<TrackerCellID> m_cellID;
    MyG4UniversalFluctuationForSi m_fluctuate;
    std::vector<std::vector<double>> m_DigitizedBins{}; // one charge bin table per layer

    mutable Gaudi::Accumulators::StatCounter<unsigned long> m_nSimHits{this, "SimTrackerHits per event"};
    mutable Gaudi::Accumulators::StatCounter<unsigned long> m_nRecoHits{this, "TrackerHits per event"};
    mutable Gaudi::Accumulators::Counter<> m_nOffSurface{this, "SimTrackerHits outside their surface"};
    mutable Gaudi::Accumulators::Counter<> m_nNoCharge{this, "SimTrackerHits without charge above threshold"};
    mutable Gaudi::Accumulators::MsgCounter<MSG::ERROR> m_unmappedLayer{this, "SimTrackerHit layer ID not found in LayerIDs"};
    mutable Gaudi::Accumulators::MsgCounter<MSG::ERROR> m_invalidTimeBinning{this,
        "Invalid setting for pixel time digitization binning. Retaining original time."};

    /* Geometry */
    void FillChargeBins();
    void PrintGeometryInfo() const;

    /* Charge digitization helpers */
    void ProduceIonisationPoints(const edm4hep::SimTrackerHit& hit, HitState& state, EventRandom& random) const;
    void ProduceSignalPoints(HitState& state) const;
    void ProduceHits(PixelHitMap& pixels, const edm4hep::SimTrackerHit& simHit, const HitState& state) const;
    void PoissonSmearer(PixelHitMap& pixels, EventRandom& random) const;
    void GainSmearer(PixelHitMap& pixels, EventRandom& random) const;
    void ApplyThreshold(PixelHitMap& pixels, EventRandom& random) const;
    void ChargeDigitizer(PixelHitMap& pixels, const HitState& state) const;

    /* Time digitization helpers */
    void TimeSmearer(PixelHitMap& pixels, const HitState& state, EventRandom& random) const;
    void TimeDigitizer(PixelHitMap& pixels) const;

    /* Reconstruction of measurement and helpers */
    bool ReconstructTrackerHit(const PixelHitMap& pixels, const HitState& state, ClusterHit& cluster) const;
    void TransformToLab(const std::uint64_t cellID, const double* xLoc, double* xLab) const;
    void FindLocalPosition(const edm4hep::SimTrackerHit& hit, double* localPosition, double* localDirection,
                           HitState& state) const;
    void TransformXYToCellID(double x, double y, int& ix, int& iy, const HitState& state) const;
    void TransformCellIDToXY(int ix, int iy, double& x, double& y, const HitState& state) const;
    int GetPixelsInaRow(const HitState& state) const;
    int GetPixelsInaColumn(const HitState& state) const;

    double randomTail(const double qmin, const double qmax, EventRandom& random) const;
    int layerMapping(int id) const;
};

#endif //MuonCVXDDigitiser_h
