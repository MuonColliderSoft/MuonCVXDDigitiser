#ifndef MuonCVXDDigitiser_h
#define MuonCVXDDigitiser_h 1

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
#include "MyG4UniversalFluctuationForSi.h"

using marlin::Processor;

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
 * @param ThresholdSmearSigma: sigma of Gaussian used in threshold smearing <br>
 * (default parameter value: 25) <br>
 * @param ChargeDigitize flag to switch on charge discretization <br>
 * (default parameter value: 1) <br>
 * @param ChargeDigitizeNumBits number of bits used to determine bins for charge discretization <br>
 * (default parameter value: 4) <br>
 * @param ChargeDigitizeBinning binning scheme for charge discretization
 * (default parameter value; 1) <br>
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
class MuonCVXDDigitiser : public Processor
{  
public:
  
    virtual Processor*  newProcessor() { return new MuonCVXDDigitiser ; }

    MuonCVXDDigitiser();

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

    int _nRun;
    int _nEvt;
    int _debug;
    int _totEntries;
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
    int _PoissonSmearing;
    int _thresholdSmearSigma;
    int _DigitizeCharge;
    int _ChargeDigitizeNumBits;
    int _ChargeDigitizeBinning;
    double _chargeMax;
    int _DigitizeTime;
    int _TimeDigitizeNumBits;
    double _timeMax;
    int _TimeDigitizeBinning;
    double _timeSmearingSigma;
    int _electronicEffects;
    int _produceFullPattern;

    MyG4UniversalFluctuationForSi *_fluctuate;

    // charge discretization
    std::vector<double> _DigitizedBins{};
    
    // geometry
    int _numberOfLayers;
    std::vector<int>   _laddersInLayer{};
#ifdef ZSEGMENTED
    std::vector<int>   _sensorsPerLadder{};
#endif
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

    // internal state
    int _currentLayer;
    int _currentLadder;
    int _numberOfSegments;
    double _currentParticleMass;
    double _currentParticleMomentum;
    double _currentPhi;
    double _eSum;
    double _segmentDepth;
    double _currentLocalPosition[3];
    double _currentEntryPoint[3];
    double _currentExitPoint[3];
    IonisationPointVec _ionisationPoints;
    SignalPointVec _signalPoints;

    /* Charge digitization helpers */
    void ProduceIonisationPoints(SimTrackerHit *hit);
    void ProduceSignalPoints();
    void ProduceHits(SimTrackerHitImplVec &simTrkVec, SimTrackerHit &simHit);
    void PoissonSmearer(SimTrackerHitImplVec &simTrkVec);
    void GainSmearer(SimTrackerHitImplVec &simTrkVec);
    void ApplyThreshold(SimTrackerHitImplVec &simTrkVec);
    void ChargeDigitizer(SimTrackerHitImplVec &simTrkVec);

    /* Time digitization helpers */
    void TimeSmearer(SimTrackerHitImplVec &simTrkVec);
    void TimeDigitizer(SimTrackerHitImplVec &simTrkVec);

    /* Reconstruction of measurement and helpers */
    TrackerHitPlaneImpl *ReconstructTrackerHit(SimTrackerHitImplVec &simTrkVec);
    void TransformToLab(const int cellID, const double *xLoc, double *xLab);
    void FindLocalPosition(SimTrackerHit *hit, double *localPosition, double *localDirection);
    void TransformXYToCellID(double x, double y, int & ix, int & iy);
    void TransformCellIDToXY(int ix, int iy, double & x, double & y);
    int GetPixelsInaRow();
    int GetPixelsInaColumn();

    void LoadGeometry();
    void PrintGeometryInfo();
    double randomTail( const double qmin, const double qmax );
};

#endif



