#ifndef MuonCVXDDigitiser_h
#define MuonCVXDDigitiser_h 1

#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "marlin/Processor.h"
#include "lcio.h"
#include "EVENT/SimTrackerHit.h"
#include "EVENT/LCIO.h"
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
 * @param CollectionName name of the input SimTrackerHit collection <br>
 * (default parameter value : "VertexBarrelCollection") <br>
 * @param OutputCollectionName name of the output TrackerHitPlane collection <br>
 * (default parameter value : "VTXTrackerHits") <br>
 * @param RelationColName name of the output LCRelation collection, reconstructed to simulated hits <br>
 * (default parameter value : "VTXTrackerHitRelations") <br>
 * @param SubDetectorName name of the sub-detector, which sets barrel or endcap and vertex,
 * inner or outer tracker <br>
 * (default parameter value : "VertexBarrel") <br>
 * @param TanLorentz tangent of the Lorentz angle along x <br>
 * (default parameter value : 0.8) <br>
 * @param TanLorentzY tangent of the Lorentz angle along y <br>
 * (default parameter value : 0) <br>
 * @param CutOnDeltaRays cut on the energy of delta-electrons (in MeV) <br>
 * (default parameter value : 0.03) <br>
 * @param DiffusionCoefficient diffusion coefficient sqrt(2D / mu / V): the diffusion sigma is this
 * times the drift distance <br>
 * (default parameter value : 0.07) <br>
 * @param PixelSizeX pixel size along u, perpendicular to the beam axis in the barrel (in mm) <br>
 * (default parameter value : 0.025) <br>
 * @param PixelSizeY pixel size along v, along the beam axis in the barrel (in mm) <br>
 * (default parameter value : 0.025) <br>
 * @param ElectronsPerKeV number of electrons produced per keV of deposited energy <br>
 * (default parameter value : 270.3) <br>
 * @param Threshold threshold on the charge collected by one pixel (in electrons) <br>
 * (default parameter value : 500) <br>
 * @param ChargeMaximum upper end of the uniform charge digitization (ChargeDigitizeBinning = 0),
 * in electrons <br>
 * (default parameter value : 15000) <br>
 * @param SegmentLength length of the segments the track is split into inside the sensor (in mm);
 * the number of segments is ceil(track length / SegmentLength) <br>
 * (default parameter value : 0.005) <br>
 * @param PoissonSmearing flag to switch on the Poisson (Gaussian above 1000 electrons) smearing
 * of the electrons collected by each pixel <br>
 * (default parameter value : 1) <br>
 * @param ThresholdSmearSigma sigma of the Gaussian smearing of the threshold of each pixel (in electrons) <br>
 * (default parameter value : 25) <br>
 * @param DigitizeCharge flag to switch on the charge digitization <br>
 * (default parameter value : 1) <br>
 * @param ChargeDigitizeNumBits number of bits of the charge digitization <br>
 * (default parameter value : 4) <br>
 * @param ChargeDigitizeBinning binning of the charge digitization: 0 = uniform, 1 = variable <br>
 * (default parameter value : 1) <br>
 * @param DigitizeTime flag to switch on the time digitization <br>
 * (default parameter value : 1) <br>
 * @param TimeDigitizeNumBits number of bits of the time digitization <br>
 * (default parameter value : 10) <br>
 * @param TimeDigitizeBinning binning of the time digitization: 0 = uniform <br>
 * (default parameter value : 0) <br>
 * @param TimeMaximum dynamic range of the time measurement (in ns) <br>
 * (default parameter value : 10) <br>
 * @param TimeSmearingSigma constant Gaussian time resolution of each pixel (in ns), used when
 * TimeSmearingModel = 1 <br>
 * (default parameter value : 0.05) <br>
 * @param TimeSmearingModel time smearing model: 0 = none, 1 = constant sigma taken from
 * TimeSmearingSigma, 2 = realistic resolution derived from the sensor thickness <br>
 * (default parameter value : 1) <br>
 * @param DoMultipleScattering Flag to enable multiple scattering of the track inside the sensor <br>
 * (default parameter value : 0) <br>
 * @param MSSliceThickness thickness (mm) of the slices the track is stepped through when
 * applying multiple scattering; the sensor is divided into a whole number of slices no
 * thicker than this <br>
 * (default parameter value : 0.005) <br>
 * @param ResimulateIonisation source of the ionisation trail: 0 = path length and deposited
 * energy from the Geant4 hit, 1 = re-simulated from the local direction and EnergyLoss <br>
 * (default parameter value : 0) <br>
 * @param SigmaLandau override for Landau sigma <br>
 * (default parameter value : -1, i.e. use the thickness-derived default) <br>
 * @param SigmaTimewalk override for timewalk sigma <br>
 * (default parameter value : -1, i.e. use the thickness-derived default) <br>
 * @param SigmaJitter override for jitter sigma <br>
 * (default parameter value : -1, i.e. use the thickness-derived default) <br>
 * @param SigmaTDC override for TDC sigma <br>
 * (default parameter value : -1, i.e. use the thickness-derived default) <br>
 * @param SigmaClock override for clock sigma <br>
 * (default parameter value : -1, i.e. use the thickness-derived default) <br>
 * @param TRise override for rise time <br>
 * (default parameter value : -1, i.e. use the thickness-derived default) <br>
 * @param ElectronicEffects flag to switch on the Gaussian smearing of the signal by the electronic noise <br>
 * (default parameter value : 1) <br>
 * @param ElectronicNoise electronic noise (in electrons) <br>
 * (default parameter value : 80) <br>
 * @param StoreFiredPixels flag to also store the fired pixels, in VBPixels, VEPixels, IBPixels,
 * IEPixels, OBPixels or OEPixels depending on SubDetectorName <br>
 * (default parameter value : 0) <br>
 * @param EnergyLoss mean energy loss along the track, around which the energy of each segment is
 * sampled (in keV/mm) <br>
 * (default parameter value : 280) <br>
 * @param MaxEnergyDelta charge added to the segments, sampled from 1/n^2, until the total is within
 * this of the Geant4 deposit (in electrons) <br>
 * (default parameter value : 100) <br>
 * @param MaxTrackLength maximum track length inside the sensor (in mm) <br>
 * (default parameter value : 10) <br>
 * @param LayerIDs IDs of the layers to digitise, in the order of the layers of the geometry;
 * hits on other layers are skipped <br>
 * (default parameter value : empty) <br>
 * <br>
 */
class MuonCVXDDigitiser : public Processor
{  
public:
  
    virtual Processor*  newProcessor() { return new MuonCVXDDigitiser ; }

    MuonCVXDDigitiser();
    ~MuonCVXDDigitiser();

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

    std::string _subDetName;
    bool isBarrel{false};
    bool isVertex{false};

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
    int _timeSmearingModel;
    int _electronicEffects;
    int _produceFullPattern;
    int _resimulateIonisation;
    int _doMultipleScattering;
    double _msSliceThickness;
    std::vector<int> _layerIDs;

    // time digitization overrides
    double _t_riseOverride;
    double _sigma_landauOverride;
    double _sigma_timewalkOverride;
    double _sigma_jitterOverride;
    double _sigma_TDCOverride;
    double _sigma_clockOverride;

    // charge discretization
    std::vector<std::vector<double>> _DigitizedBins{};  // one bin table per layer
    
    // geometry
    int _numberOfLayers;
    std::vector<float> _layerThickness{};
    std::vector<float> _layerHalfThickness{};
    const dd4hep::rec::SurfaceMap* _map ;
    // Half-lengths along u and v (mm) of the sensitive area of each sensor, by cell ID.
    // The pixel matrix covers the sensor the hit is on: the ladder and petal dimensions
    // of the layer do not match all the sensors (e.g. tracker endcap modules).
    std::unordered_map<unsigned long, std::pair<double, double>> _sensorHalfLengths{};

    // internal state
    double _sensorHalfLengthU{0};
    double _sensorHalfLengthV{0};
    int _currentLayer;
    int _currentLadder;
    int _numberOfSegments;
    double _currentParticleMass;
    double _currentParticleMomentum;
    double _eSum;
    std::vector<std::array<double, 3>> _msTrajectory;
    std::vector<double> _msSegmentLengths;
    std::vector<double> _msCumulativePath;
    double _msTotalScatteringAngle;
    IonisationPointVec _ionisationPoints;
    SignalPointVec _signalPoints;
    MyG4UniversalFluctuationForSi *_fluctuate;

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
    int layerMapping( int id);
};

#endif



