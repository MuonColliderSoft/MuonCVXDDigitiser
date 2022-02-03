#ifndef MuonCVXDRealDigitiser_h
#define MuonCVXDRealDigitiser_h 1

#include <string>
#include <vector>

#include "marlin/Processor.h"
#include "lcio.h"
#include "EVENT/SimTrackerHit.h"
#include "EVENT/LCIO.h"
#include "IMPL/TrackerHitImpl.h"
#include <IMPL/TrackerHitPlaneImpl.h>
#include "IMPL/SimTrackerHitImpl.h"
#include <IMPL/LCCollectionVec.h>
#include "DDRec/Surface.h"
#include "DDRec/SurfaceManager.h"

#include "AIDA/IHistogram1D.h"

using marlin::Processor;
using AIDA::IHistogram1D;

struct IonisationPoint
{
    double x;
    double y;
    double z;
    double eloss;
};

struct SignalPoint
{
    double x;
    double y;
    double sigmaX;
    double sigmaY;
    double charge;
};

typedef std::vector<SimTrackerHitImpl*> SimTrackerHitImplVec;
typedef std::vector<IonisationPoint> IonisationPointVec;
typedef std::vector<SignalPoint> SignalPointVec;

/**  Digitizer for Simulated Hits in the Vertex Detector. <br>
 * Digitization follows the procedure adopted in the CMS software package. 
 * See https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePixelDigitization
 * 
 * @param CollectionName name of input SimTrackerHit collection <br>
 * (default parameter value : "VXDCollection")
 * @param OutputCollectionName name of output TrackerHitsPlane collection <br>
 * (default parameter value : "VTXTrackerHits")
 * @param RelationColName name of the LCRelation <br>
 * (default parameter value : "VTXTrackerHitRelations")
 * @param SubDetectorName name of the detector <br>
 * (default parameter value : "VertexBarrel")
 * @param TanLorentz tangent of the Lorentz angle <br>
 * (default parameter value : 0.8) <br>
 * @param TanLorentzY tangent of the Lorentz angle along Y <br>
 * (default parameter value : 0) <br>
 * @param CutOnDeltaRays cut on the energy of delta-electrons (in MeV) <br>
 * (default parameter value : 0.03) <br>
 * @param Diffusion diffusion coefficient for the nominal active layer thickness (in mm) <br>
 * (default parameter value : 0.002) <br>
 * @param PixelSizeX pixel size along direction perpendicular to beam axis (in mm) <br>
 * (default value : 0.025) <br>
 * @param PixelSizeY pixel size along beam axis (in mm) <br>
 * (default value : 0.025) <br>
 * @param ElectronsPerMeV number of electrons produced per MeV of deposited energy <br>
 * (default parameter value : 270.3) <br>
 * @param Threshold threshold on charge deposited on one pixel (in electons) <br>
 * (default parameter value : 200.0) <br>
 * @param SegmentLength segment length along track path which is used to subdivide track into segments (in mm).
 * The number of track subsegments is calculated as int(TrackLengthWithinActiveLayer/SegmentLength)+1 <br>
 * (default parameter value : 0.005) <br>
 * @param WidthOfCluster defines width in Gaussian sigmas to perform charge integration for 
 * a given pixel <br>
 * (default parameter value : 3.0) <br>
 * @param PoissonSmearing flag to switch on gaussian smearing of electrons collected on pixels <br>
 * (default parameter value : 1) <br>
 * @param ElectronicEffects flag to switch on gaussian smearing of signal (electronic noise) <br>
 * (default parameter value : 1) <br>
 * @param ElectronicNoise electronic noise in electrons <br>
 * (default parameter value : 100) <br>
 * @param StoreFiredPixels flag to store also the fired pixels (collection names: "VTXPixels") <br>
 * (default parameter value : 0) <br>
 * @param EnergyLoss Energy loss in keV/mm <br>
 * (default parameter value : 280.0) <br>
 * @param MaxEnergyDelta max delta in energy hit (difference from the hit simulated charged and the comuputed ones) in electrons <br>
 * (default parameter value : 100) <br>
  * @param MaxTrackLength Maximum values for track path length inside the ladder (in mm)", <br>
 * (default parameter value : 10) <br> 
 * <br>
 */
class MuonCVXDRealDigitiser : public Processor
{  
public:
  
    virtual Processor*  newProcessor() { return new MuonCVXDRealDigitiser ; }

    MuonCVXDRealDigitiser();

    /** Called at the begin of the job before anything is read.
    * Use to initialize the processor, e.g. book histograms.
    */
    virtual void init();

    /** Called for every run.
    */
    virtual void processRunHeader( LCRunHeader* run );

    /** Called for every event - the working horse.
    */
    virtual void processEvent( LCEvent * evt ); 

    virtual void check( LCEvent * evt );

    /** Called after data processing for clean up.
    */
    virtual void end();

protected:

    void PrintGeometryInfo();

    int _nRun;
    int _nEvt;
    int _debug;
    int _totEntries;
    int _barrelID;
    std::string _subDetName;

    // input/output collections
    std::string _colName;
    std::string _outputCollectionName;
    std::string _colVTXRelation;

    // processor parameters
    double _tanLorentzAngleX;
    double _tanLorentzAngleY;
    double _cutOnDeltaRays;
    double _diffusionCoefficient;
    double _pixelSizeX;
    double _pixelSizeY;
    double _electronsPerKeV;
    double _segmentLength;
    double _threshold;
    double _electronicNoise;
    double _energyLoss;
    double _deltaEne;
  	double _maxTrkLen;
    float _window_size;
    float _fe_slope;
    int _PoissonSmearing;
    int _electronicEffects;
    int _produceFullPattern;
    int sensor_type;

    // geometry
    int _numberOfLayers;
    std::vector<int>   _laddersInLayer{};
    std::vector<int>   _sensorsPerLadder{};
    std::vector<float> _layerRadius{};
    std::vector<float> _layerThickness{};
    std::vector<float> _layerHalfThickness{};
    std::vector<float> _layerLadderLength{};
    std::vector<float> _layerLadderHalfWidth{};
    std::vector<float> _layerPhiOffset{};
    std::vector<float> _layerActiveSiOffset{};
    std::vector<float> _layerHalfPhi{};
    std::vector<float> _layerLadderWidth{};
    const dd4hep::rec::SurfaceMap* _map ;

    bool create_stats;
    IHistogram1D* signal_dHisto;
    IHistogram1D* bib_dHisto;
    IHistogram1D* signal_cSizeHisto;
    IHistogram1D* signal_xSizeHisto;
    IHistogram1D* signal_ySizeHisto;
    IHistogram1D* signal_zSizeHisto;
    IHistogram1D* signal_eDepHisto;
    IHistogram1D* bib_cSizeHisto;
    IHistogram1D* bib_xSizeHisto;
    IHistogram1D* bib_ySizeHisto;
    IHistogram1D* bib_zSizeHisto;
    IHistogram1D* bib_eDepHisto;

};

#endif //MuonCVXDRealDigitiser_h



