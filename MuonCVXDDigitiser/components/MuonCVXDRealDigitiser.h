#ifndef MuonCVXDRealDigitiser_h
#define MuonCVXDRealDigitiser_h 1

#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include "Gaudi/Accumulators.h"
#include "Gaudi/Accumulators/StaticRootHistogram.h"
#include "Gaudi/Property.h"
#include "k4FWCore/Transformer.h"
#include "k4Interface/IGeoSvc.h"
#include "k4Interface/IUniqueIDGenSvc.h"

#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/TrackerHitPlaneCollection.h"
#include "edm4hep/TrackerHitSimTrackerHitLinkCollection.h"

#include "G4UniversalFluctuation.h"
#include "LayerGeometry.h"
#include "TrackerCellID.h"

/** Digitizer for Simulated Hits in the barrel of the Vertex Detector, with a time-resolved
 * model of the readout chip. <br>
 * Digitization follows the procedure adopted in the CMS software package.
 * See https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePixelDigitization
 *
 * The SimTrackerHits of every ladder are fed in time order to a sensor model through a
 * sliding time window of WindowSize ns. Pixels are clustered on each sensor of the ladder
 * as the clock advances and a TrackerHitPlane is produced for each cluster.
 *
 * Inputs:
 * - CollectionName: SimTrackerHits of the sub-detector
 * - EventHeader: used to seed the random numbers for each event
 *
 * Outputs:
 * - OutputCollectionName: reconstructed TrackerHitPlanes
 * - RelationColName: links from each reconstructed hit to all the SimTrackerHits contributing to its cluster (weight 1)
 *
 * With CreateStats, cluster size, cluster energy and hit offset histograms are filled separately
 * for clusters with at least one non-overlay SimTrackerHit (signal) and for the others (BIB).
 *
 * All the parameters are documented in the property declarations below.
 */
struct MuonCVXDRealDigitiser final
    : k4FWCore::MultiTransformer<std::tuple<edm4hep::TrackerHitPlaneCollection,
                                            edm4hep::TrackerHitSimTrackerHitLinkCollection>(
          const edm4hep::SimTrackerHitCollection&, const edm4hep::EventHeaderCollection&)>
{
    MuonCVXDRealDigitiser(const std::string& name, ISvcLocator* svcLoc);

    StatusCode initialize() override;

    std::tuple<edm4hep::TrackerHitPlaneCollection, edm4hep::TrackerHitSimTrackerHitLinkCollection>
    operator()(const edm4hep::SimTrackerHitCollection& simTrackerHits,
               const edm4hep::EventHeaderCollection& headers) const override;

private:
    using Histogram = Gaudi::Accumulators::StaticRootHistogram<1>;

    Gaudi::Property<std::string> m_subDetName{this, "SubDetectorName", "VertexBarrel", "Name of Vertex detector"};
    Gaudi::Property<std::vector<int>> m_layerIDs{this, "LayerIDs", {},
        "ID of layers of subdetector, in the order of the geometry layers. Empty = 0, 1, 2, ..."};
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
    Gaudi::Property<double> m_pixelSizeX{this, "PixelSizeX", 0.025, "Pixel Size X"};
    Gaudi::Property<double> m_pixelSizeY{this, "PixelSizeY", 0.025, "Pixel Size Y"};
    Gaudi::Property<double> m_electronsPerKeV{this, "ElectronsPerKeV", 270.3, "Electrons per keV"};
    Gaudi::Property<double> m_threshold{this, "Threshold", 200., "Cell Threshold in electrons"};
    Gaudi::Property<double> m_segmentLength{this, "SegmentLength", 0.005, "Segment Length in mm"};
    Gaudi::Property<double> m_energyLoss{this, "EnergyLoss", 280.0, "Energy Loss keV/mm"};
    Gaudi::Property<double> m_deltaEne{this, "MaxEnergyDelta", 100.0,
        "Max delta in energy between G4 prediction and random sampling for each hit in electrons"};
    Gaudi::Property<double> m_maxTrkLen{this, "MaxTrackLength", 10.0, "Maximum values for track length (in mm)"};
    Gaudi::Property<float> m_window_size{this, "WindowSize", 25., "Window size (in nsec)"};
    Gaudi::Property<float> m_fe_slope{this, "RD53Aslope", 0.1, "ADC slope for chip RD53A"};
    Gaudi::Property<int> m_sensor_type{this, "SensorType", 1, "Sensor model to be used (0 : ChipRD53A, 1 : Trivial)"};
    Gaudi::Property<int> m_zSegmented{this, "ZSegmented", -1,
        "Sensor segmentation along z, barrel layers only: -1 = auto (on for the vertex barrel), 0 = off, 1 = on."};
    Gaudi::Property<bool> m_create_stats{this, "CreateStats", false, "Fill cluster statistics histograms"};

    SmartIF<IGeoSvc> m_geoSvc;
    SmartIF<IUniqueIDGenSvc> m_uidSvc;

    // Set in initialize() and read-only afterwards
    LayerGeometry m_geo;
    std::unique_ptr<TrackerCellID> m_cellID;
    G4UniversalFluctuation m_fluctuate;

    mutable Gaudi::Accumulators::StatCounter<unsigned long> m_nSimHits{this, "SimTrackerHits per event"};
    mutable Gaudi::Accumulators::StatCounter<unsigned long> m_nRecoHits{this, "TrackerHits per event"};
    mutable Gaudi::Accumulators::StatCounter<unsigned long> m_nSimPerReco{this, "SimTrackerHits per TrackerHit"};
    mutable Gaudi::Accumulators::Counter<> m_nOffSurface{this, "SimTrackerHits outside their surface"};
    mutable Gaudi::Accumulators::MsgCounter<MSG::ERROR> m_pixelNumberError{this, "Pixel number error for ladder"};
    mutable Gaudi::Accumulators::MsgCounter<MSG::ERROR> m_segmentNumberError{this, "Segment number error for ladder"};
    mutable Gaudi::Accumulators::MsgCounter<MSG::ERROR> m_missingSurface{this, "No surface for the cell ID of a cluster"};

    std::unique_ptr<Histogram> m_signal_dHisto;
    std::unique_ptr<Histogram> m_bib_dHisto;
    std::unique_ptr<Histogram> m_signal_cSizeHisto;
    std::unique_ptr<Histogram> m_signal_xSizeHisto;
    std::unique_ptr<Histogram> m_signal_ySizeHisto;
    std::unique_ptr<Histogram> m_signal_zSizeHisto;
    std::unique_ptr<Histogram> m_signal_eDepHisto;
    std::unique_ptr<Histogram> m_bib_cSizeHisto;
    std::unique_ptr<Histogram> m_bib_xSizeHisto;
    std::unique_ptr<Histogram> m_bib_ySizeHisto;
    std::unique_ptr<Histogram> m_bib_zSizeHisto;
    std::unique_ptr<Histogram> m_bib_eDepHisto;

    void PrintGeometryInfo() const;
};

#endif //MuonCVXDRealDigitiser_h
