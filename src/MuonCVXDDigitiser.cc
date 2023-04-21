#include "MuonCVXDDigitiser.h"
#include <iostream>
#include <algorithm>

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>

#include <UTIL/CellIDEncoder.h>
#include <UTIL/CellIDDecoder.h>
#include "UTIL/LCTrackerConf.h"

#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCRelationImpl.h>

#include "DD4hep/Detector.h"
#include "DDRec/DetectorData.h"
#include "DD4hep/DD4hepUnits.h"

#include "gsl/gsl_sf_erf.h"
#include "gsl/gsl_math.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandPoisson.h" 
#include "CLHEP/Random/RandFlat.h" 
    
// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

using CLHEP::RandGauss;
using CLHEP::RandPoisson;
using CLHEP::RandFlat;

using dd4hep::Detector;
using dd4hep::DetElement;
using dd4hep::rec::ZPlanarData;
using dd4hep::rec::SurfaceManager;
using dd4hep::rec::SurfaceMap;
using dd4hep::rec::ISurface;
using dd4hep::rec::Vector2D;
using dd4hep::rec::Vector3D;

MuonCVXDDigitiser aMuonCVXDDigitiser ;

MuonCVXDDigitiser::MuonCVXDDigitiser() :
    Processor("MuonCVXDDigitiser"),
    _nRun(0),
    _nEvt(0),
    _totEntries(0),
    _currentLayer(0)
{
    _description = "MuonCVXDDigitiser should create VTX TrackerHits from SimTrackerHits";

    registerInputCollection(LCIO::SIMTRACKERHIT,
                           "CollectionName", 
                           "Name of the SimTrackerHit collection",
                           _colName,
                           std::string("VertexBarrelCollection"));

    registerOutputCollection(LCIO::TRACKERHITPLANE,
                            "OutputCollectionName", 
                            "Name of the output TrackerHit collection",
                            _outputCollectionName,
                            std::string("VTXTrackerHits"));

    registerOutputCollection(LCIO::LCRELATION,
                            "RelationColName", 
                            "Name of the output VTX trackerhit relation collection",
                            _colVTXRelation,
                            std::string("VTXTrackerHitRelations"));

    registerProcessorParameter("SubDetectorName", 
                               "Name of Vertex detector",
                               _subDetName,
                               std::string("VertexBarrel"));

    registerProcessorParameter("TanLorentz",
                               "Tangent of Lorentz Angle",
                               _tanLorentzAngleX,
                               (double)0.8);

    registerProcessorParameter("TanLorentzY",
                               "Tangent of Lorentz Angle along Y",
                               _tanLorentzAngleY,
                               (double)0);

    registerProcessorParameter("CutOnDeltaRays",
                               "Cut on delta-ray energy (MeV)",
                               _cutOnDeltaRays,
                               (double)0.030);

    // For diffusion-coeffieint calculation, see e.g. https://www.slac.stanford.edu/econf/C060717/papers/L008.PDF
    // or directly Eq. 13 of https://cds.cern.ch/record/2161627/files/ieee-tns-07272141.pdf
    // diffusionCoefficient = sqrt(2*D / mu / V), where
    //  - D = 12 cm^2/s // diffusion constant
    //  - mu = 450 cm^2/s/V // mobility
    //  - V = 10-30 V // expected depletion voltage
    //  => _diffusionCoefficient = 0.04-0.07
    registerProcessorParameter("DiffusionCoefficient",
                               "Diffusion coefficient, sqrt(D / mu / V).",
                               _diffusionCoefficient,
                               (double)0.07);

    registerProcessorParameter("PixelSizeX",
                               "Pixel Size X",
                               _pixelSizeX,
                               (double)0.025);

    registerProcessorParameter("PixelSizeY",
                               "Pixel Size Y",
                               _pixelSizeY,
                               (double)0.025);

    registerProcessorParameter("ElectronsPerKeV",
                               "Electrons per keV",
                               _electronsPerKeV,
                               (double)270.3);

    registerProcessorParameter("Threshold",
                               "Cell Threshold in electrons",
                               _threshold,
                               500.);

    registerProcessorParameter("ChargeMaximum",
                               "Cell dynamic range in electrons",
                               _chargeMax,
                               15000.);

    registerProcessorParameter("SegmentLength",
                               "Segment Length in mm",
                               _segmentLength,
                               double(0.005));

    registerProcessorParameter("PoissonSmearing",
                               "Apply Poisson smearing of electrons collected on pixels",
                               _PoissonSmearing,
                               1);
    
    registerProcessorParameter("ThresholdSmearSigma",
                               "sigma of Gaussian used in threshold smearing, in electrons",
                               _thresholdSmearSigma, 
                               25);
    
    registerProcessorParameter("DigitizeCharge",
                               "Flag to enable Digitization of the charge collected on pixels",
                               _DigitizeCharge,
                               1);
    
    registerProcessorParameter("ChargeDigitizeNumBits",
                               "Number of bits used to determine bins for charge discretization",
                               _ChargeDigitizeNumBits,
                               4);

    registerProcessorParameter("ChargeDigitizeBinning",
                               "Binning scheme used for charge discretization",
                               _ChargeDigitizeBinning,
                               1);

    registerProcessorParameter("DigitizeTime",
                                "Flag to enable digitization of timing information.",
                                _DigitizeTime, 
                                1);

    registerProcessorParameter("TimeDigitizeNumBits",
                               "Number of bits used to determine bins for time discretization",
                               _TimeDigitizeNumBits,
                               10);

    registerProcessorParameter("TimeDigitizeBinning",
                               "Binning scheme used for time discretization",
                               _TimeDigitizeBinning,
                               0);

    registerProcessorParameter("TimeMaximum",
                               "Cell dynamic range for timing measurement [ns]",
                               _timeMax,
                               10.000);

    registerProcessorParameter("TimeSmearingSigma",
                                "Effective intrinsic time measurement resolution effects (ns).",
                                _timeSmearingSigma, 
                                0.05);

    registerProcessorParameter("ElectronicEffects",
                               "Apply Electronic Effects",
                               _electronicEffects,
                               int(1));

    registerProcessorParameter("ElectronicNoise",
                               "electronic noise in electrons",
                               _electronicNoise,
                               80.);

    registerProcessorParameter("StoreFiredPixels",
                               "Store fired pixels",
                               _produceFullPattern,
                               int(0));

    registerProcessorParameter("EnergyLoss",
                               "Energy Loss keV/mm",
                               _energyLoss,
                               double(280.0));
                               
    registerProcessorParameter("MaxEnergyDelta",
                               "Max delta in energy between G4 prediction and random sampling for each hit in electrons",
                               _deltaEne,
                               100.0);                               

    registerProcessorParameter("MaxTrackLength",
                               "Maximum values for track length (in mm)",
                               _maxTrkLen,
                               10.0); 

}



void MuonCVXDDigitiser::init()
{ 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;

    printParameters() ;

    _nRun = 0 ;
    _nEvt = 0 ;
    _totEntries = 0;
    _fluctuate = new MyG4UniversalFluctuationForSi();
}


void MuonCVXDDigitiser::processRunHeader(LCRunHeader* run)
{ 
    _nRun++ ;

    LoadGeometry() ;
}


void MuonCVXDDigitiser::LoadGeometry()
{
    Detector& theDetector = Detector::getInstance();
    DetElement vxBarrel = theDetector.detector(_subDetName);              // TODO check missing barrel
    ZPlanarData&  zPlanarData = *vxBarrel.extension<ZPlanarData>();       // TODO check missing extension
    std::vector<ZPlanarData::LayerLayout> vx_layers = zPlanarData.layers;
    _numberOfLayers  = vx_layers.size();

    SurfaceManager& surfMan = *theDetector.extension<SurfaceManager>();
    _map = surfMan.map( vxBarrel.name() ) ;
    if( ! _map ) 
    {
      std::stringstream err  ; err << " Could not find surface map for detector: "
                                 << _subDetName << " in SurfaceManager " ;
      throw Exception( err.str() ) ;
    }

    _laddersInLayer.resize(_numberOfLayers);
#ifdef ZSEGMENTED
    _sensorsPerLadder.resize(_numberOfLayers);
#endif
    _layerHalfPhi.resize(_numberOfLayers);
    _layerHalfThickness.resize(_numberOfLayers);
    _layerThickness.resize(_numberOfLayers);
    _layerRadius.resize(_numberOfLayers);
    _layerLadderLength.resize(_numberOfLayers);
    _layerLadderWidth.resize(_numberOfLayers);
    _layerLadderHalfWidth.resize(_numberOfLayers);
    _layerActiveSiOffset.resize(_numberOfLayers);
    _layerPhiOffset.resize(_numberOfLayers);

    int curr_layer = 0;
    for(ZPlanarData::LayerLayout z_layout : vx_layers)
    {
        // ALE: Geometry is in cm, convert all lenght in mm
        _laddersInLayer[curr_layer] = z_layout.ladderNumber;

        _layerHalfPhi[curr_layer] = M_PI / ((double)_laddersInLayer[curr_layer]) ;

        _layerThickness[curr_layer] = z_layout.thicknessSensitive * dd4hep::cm / dd4hep::mm ;

        _layerHalfThickness[curr_layer] = 0.5 * _layerThickness[curr_layer];

        _layerRadius[curr_layer] = z_layout.distanceSensitive * dd4hep::cm / dd4hep::mm  + _layerHalfThickness[curr_layer];

#ifdef ZSEGMENTED
        _sensorsPerLadder[curr_layer] = z_layout.sensorsPerLadder;

        _layerLadderLength[curr_layer] = z_layout.lengthSensor * z_layout.sensorsPerLadder * dd4hep::cm / dd4hep::mm ;
#else
        _layerLadderLength[curr_layer] = z_layout.lengthSensor * dd4hep::cm / dd4hep::mm ;
#endif
        _layerLadderWidth[curr_layer] = z_layout.widthSensitive * dd4hep::cm / dd4hep::mm ;

        _layerLadderHalfWidth[curr_layer] = _layerLadderWidth[curr_layer] / 2.;

        _layerActiveSiOffset[curr_layer] = - z_layout.offsetSensitive * dd4hep::cm / dd4hep::mm ;

        _layerPhiOffset[curr_layer] = z_layout.phi0;

        curr_layer++;
    }

    PrintGeometryInfo();

    // Bins for charge discretization
    // FIXME: Will move to assign more dynamically 
    if (_ChargeDigitizeNumBits == 3) _DigitizedBins = {500, 786, 1100, 1451, 1854, 2390, 3326, 31973};
    if (_ChargeDigitizeNumBits == 4) _DigitizedBins = {500, 639, 769, 910, 1057, 1213, 1379, 1559, 1743, 1945, 2193, 2484, 2849, 3427, 4675, 29756};
    if (_ChargeDigitizeNumBits == 5) _DigitizedBins = {500, 573, 633, 698, 757, 821, 890, 963, 1032, 1104, 1179, 1260, 1337, 1421, 1505, 1600, 1685, 1777, 1875, 1982, 2097, 2220, 2352, 2511, 2679, 2866, 3107, 3429, 3880, 4618, 6287, 16039};
    if (_ChargeDigitizeNumBits == 6) _DigitizedBins = {500, 542, 572, 601, 629, 661, 692, 721, 750, 779, 812, 842, 877, 913, 946, 981, 1016, 1051, 1087, 1121, 1161, 1196, 1237, 1275, 1313, 1350, 1391, 1431, 1468, 1514, 1560, 1606, 1646, 1687, 1733, 1777, 1821, 1872, 1920, 1976, 2036, 2091, 2145, 2213, 2272, 2337, 2411, 2488, 2573, 2651, 2739, 2834, 2938, 3053, 3194, 3356, 3532, 3764, 4034, 4379, 4907, 5698, 6957, 9636};
    if (_ChargeDigitizeNumBits == 8) _DigitizedBins = {500, 511, 523, 533, 542, 550, 556, 564, 570, 577, 585, 592, 598, 603, 610, 617, 624, 630, 638, 646, 654, 661, 668, 676, 684, 691, 699, 705, 712, 719, 724, 731, 738, 745, 752, 760, 767, 772, 780, 787, 795, 802, 810, 818, 826, 832, 839, 847, 856, 865, 874, 881, 889, 898, 906, 916, 924, 930, 938, 945, 955, 965, 971, 978, 986, 995, 1004, 1012, 1019, 1027, 1036, 1044, 1053, 1062, 1071, 1079, 1088, 1096, 1104, 1112, 1121, 1131, 1139, 1149, 1158, 1168, 1175, 1184, 1193, 1203, 1211, 1221, 1233, 1241, 1249, 1259, 1268, 1277, 1286, 1294, 1303, 1313, 1321, 1330, 1338, 1348, 1357, 1368, 1378, 1387, 1395, 1406, 1417, 1426, 1434, 1445, 1452, 1460, 1470, 1480, 1492, 1503, 1514, 1525, 1536, 1550, 1560, 1570, 1580, 1592, 1604, 1614, 1623, 1634, 1644, 1653, 1662, 1673, 1684, 1695, 1707, 1717, 1727, 1737, 1747, 1759, 1769, 1780, 1790, 1800, 1812, 1823, 1835, 1846, 1860, 1873, 1885, 1897, 1907, 1918, 1931, 1943, 1958, 1971, 1987, 2000, 2014, 2026, 2041, 2056, 2068, 2080, 2095, 2108, 2119, 2131, 2147, 2162, 2180, 2195, 2213, 2224, 2238, 2256, 2269, 2284, 2300, 2314, 2332, 2351, 2366, 2383, 2401, 2421, 2440, 2458, 2475, 2496, 2519, 2538, 2559, 2581, 2601, 2618, 2636, 2658, 2681, 2703, 2722, 2742, 2767, 2791, 2811, 2836, 2857, 2884, 2913, 2938, 2967, 2995, 3023, 3052, 3086, 3119, 3153, 3188, 3221, 3270, 3304, 3342, 3390, 3428, 3473, 3515, 3556, 3611, 3691, 3742, 3801, 3857, 3928, 3999, 4069, 4141, 4220, 4325, 4417, 4518, 4655, 4789, 4965, 5141, 5359, 5548, 5770, 6017, 6311, 6584, 7024, 7492, 8060, 8740, 9738, 11450, 14878, 23973};

} 



void MuonCVXDDigitiser::processEvent(LCEvent * evt)
{ 

    //SP. few TODO items:
    // - include noisy pixels (calculate rate from gaussian with unit sigma integral x > _electronicNoise / _threshold )
    // - change logic in creating pixels from all SimTrkHits, then cluster them (incl. timing info)
    // - include threshold dispersion effects
    // - add digi parametrization for time measurement
    // - change position determination of cluster to analog cluster (w-avg of corner hits)

    if ((_nEvt == 0) and (!_map)){
        LoadGeometry() ;
    }

    LCCollection * STHcol = nullptr;
    try
    {
        STHcol = evt->getCollection(_colName);
    }
    catch( lcio::DataNotAvailableException &ex )
    {
        streamlog_out(WARNING) << _colName << " collection not available" << std::endl;
        STHcol = nullptr;
    }

    if( STHcol != nullptr )
    {
        LCCollectionVec *THcol = new LCCollectionVec(LCIO::TRACKERHITPLANE);
        CellIDEncoder<TrackerHitPlaneImpl> cellid_encoder( lcio::LCTrackerCellID::encoding_string(), THcol ) ;

        CellIDDecoder<SimTrackerHit> cellid_decoder( STHcol) ;

        LCCollectionVec* relCol = new LCCollectionVec(LCIO::LCRELATION);
       // to store the weights
       LCFlagImpl lcFlag(0) ;
       lcFlag.setBit( LCIO::LCREL_WEIGHTED ) ;
       relCol->setFlag( lcFlag.getFlag()  ) ;

        LCCollectionVec *STHLocCol = nullptr;
        if (_produceFullPattern != 0)
        {
            STHLocCol = new LCCollectionVec(LCIO::SIMTRACKERHIT);
            CellIDEncoder<SimTrackerHitImpl> cellid_encoder_fired( lcio::LCTrackerCellID::encoding_string(), STHLocCol ) ;
        }

        int nSimHits = STHcol->getNumberOfElements();

        streamlog_out( DEBUG9 ) << "Processing collection " << _colName  << " with " <<  nSimHits  << " hits ... " << std::endl ;

        for (int i=0; i < nSimHits; ++i)
        {
            SimTrackerHit * simTrkHit = 
                dynamic_cast<SimTrackerHit*>(STHcol->getElementAt(i));

            // use CellID to set layer and ladder numbers
            _currentLayer  = cellid_decoder( simTrkHit )["layer"];
            _currentLadder = cellid_decoder( simTrkHit )["module"];
            streamlog_out( DEBUG7 ) << "Processing simHit #" << i << ", from layer=" << _currentLayer << ", module=" << _currentLadder << std::endl;
            streamlog_out (DEBUG6) << "- EDep = " << simTrkHit->getEDep() *dd4hep::GeV / dd4hep::keV << " keV, path length = " << simTrkHit->getPathLength() * 1000. << " um" << std::endl;
            float mcp_r = std::sqrt(simTrkHit->getPosition()[0]*simTrkHit->getPosition()[0]+simTrkHit->getPosition()[1]*simTrkHit->getPosition()[1]);
            float mcp_phi = std::atan(simTrkHit->getPosition()[1]/simTrkHit->getPosition()[0]);
            float mcp_theta = simTrkHit->getPosition()[2] == 0 ? 3.1416/2 : std::atan(mcp_r/simTrkHit->getPosition()[2]);
            streamlog_out (DEBUG6) << "- Position (mm) x,y,z,t = " << simTrkHit->getPosition()[0] << ", " << simTrkHit->getPosition()[1] << ", " << simTrkHit->getPosition()[2] << ", " << simTrkHit->getTime() << std::endl;
            streamlog_out (DEBUG6) << "- Position r(mm),phi,theta = " << mcp_r << ", " << mcp_phi << ", " << mcp_theta << std::endl;
            streamlog_out (DEBUG6) << "- MC particle pdg = ";
            EVENT::MCParticle *mcp = simTrkHit->getMCParticle();
            if (mcp) {
                streamlog_out (DEBUG6) << simTrkHit->getMCParticle()->getPDG();
            } else {
                streamlog_out (DEBUG6) << " N.A.";
            }
            streamlog_out (DEBUG6) << std::endl;
            streamlog_out (DEBUG6) << "- MC particle p (GeV) = " << std::sqrt(simTrkHit->getMomentum()[0]*simTrkHit->getMomentum()[0]+simTrkHit->getMomentum()[1]*simTrkHit->getMomentum()[1]+simTrkHit->getMomentum()[2]*simTrkHit->getMomentum()[2]) << std::endl;

            streamlog_out (DEBUG6) << "- isSecondary = " << simTrkHit->isProducedBySecondary() << ", isOverlay = " << simTrkHit->isOverlay() << std::endl;
            streamlog_out (DEBUG6) << "- Quality = " << simTrkHit->getQuality() << std::endl;

            ProduceIonisationPoints( simTrkHit );       
            if (_currentLayer == -1)
              continue;

            ProduceSignalPoints();

            SimTrackerHitImplVec simTrkHitVec;

            ProduceHits(simTrkHitVec, *simTrkHit);

            if (_PoissonSmearing != 0) PoissonSmearer(simTrkHitVec);

            if (_electronicEffects != 0) GainSmearer(simTrkHitVec);

	        ApplyThreshold(simTrkHitVec);

	        if (_DigitizeCharge != 0) ChargeDigitizer(simTrkHitVec);

            if (_timeSmearingSigma > 0) TimeSmearer(simTrkHitVec);

            if (_DigitizeTime != 0) TimeDigitizer(simTrkHitVec);
	    

            //**************************************************************************
            // Create reconstructed cluster object (TrackerHitImpl)
            //**************************************************************************
            TrackerHitPlaneImpl *recoHit = ReconstructTrackerHit(simTrkHitVec);
            if (recoHit == nullptr)
            {
                streamlog_out(DEBUG) << "Skip hit" << std::endl;
                continue;
            }       
            // hit's layer/ladder position does not change
            const int cellid0 = simTrkHit->getCellID0();
            const int cellid1 = simTrkHit->getCellID1();

            recoHit->setCellID0( cellid0 );
            recoHit->setCellID1( cellid1 );

            double xLab[3];
            TransformToLab( cellid0, recoHit->getPosition(), xLab);
            recoHit->setPosition( xLab );

            SurfaceMap::const_iterator sI = _map->find( cellid0 ) ;
            const dd4hep::rec::ISurface* surf = sI->second ;

            dd4hep::rec::Vector3D u = surf->u() ;
            dd4hep::rec::Vector3D v = surf->v() ;
            
            double localPos[3];
            double localIdx[3];
            double localDir[3];
            FindLocalPosition(simTrkHit, localPos, localDir);
            localIdx[0] = localPos[0] / _pixelSizeX;
            localIdx[1] = localPos[1] / _pixelSizeY;
            float incidentPhi = std::atan(localDir[0] / localDir[2]);
            float incidentTheta = std::atan(localDir[1] / localDir[2]);
            
            float u_direction[2] ;
            //TODO HACK: Store incidence angle of particle instead!
            u_direction[0] = u.theta();
            u_direction[1] = u.phi();

            float v_direction[2] ;
            v_direction[0] = v.theta();
            v_direction[1] = v.phi();

            recoHit->setU( u_direction ) ;
            recoHit->setV( v_direction ) ;
            
            //**************************************************************************
            // Set Relation to SimTrackerHit
            //**************************************************************************    

            LCRelationImpl* rel = new LCRelationImpl;
            rel->setFrom (recoHit);
            rel->setTo (simTrkHit);
            rel->setWeight( 1.0 );
            relCol->addElement(rel);


            streamlog_out (DEBUG7) << "Reconstructed pixel cluster:" << std::endl;
            streamlog_out (DEBUG7) << "- local position (x,y) = " << localPos[0] << "(Idx: " << localIdx[0] << "), " << localPos[1] << "(Idy: " << localIdx[1] << ")" << std::endl;
            streamlog_out (DEBUG7) << "- global position (x,y,z, t) = " << recoHit->getPosition()[0] << ", " << recoHit->getPosition()[1] << ", " << recoHit->getPosition()[2] << ", " << recoHit->getTime() << std::endl;
            streamlog_out (DEBUG7) << "- charge = " << recoHit->getEDep() << "(True: " << simTrkHit->getEDep() << ")"  << std::endl;
            streamlog_out (DEBUG7) << "- incidence angles: theta = " << incidentTheta << ", phi = " << incidentPhi << std::endl;


            if (_produceFullPattern != 0)
            {
              // Store all the fired points
              for (int iS = 0; iS < (int)simTrkHitVec.size(); ++iS)
              {
                SimTrackerHitImpl *sth = simTrkHitVec[iS];
                float charge = sth->getEDep();
                //store hits that are above threshold. In case of _ChargeDiscretization, just check for a small non-zero value                
                if ( ((_DigitizeCharge > 0) and (charge >1.0)) or
                    (charge > _threshold) )
                {
                   SimTrackerHitImpl *newsth = new SimTrackerHitImpl();
                   // hit's layer/ladder position is the same for all fired points 
                   newsth->setCellID0( cellid0 );
                   newsth->setCellID1( cellid1 );
                   //Store local position in units of pixels instead
                   const double *sLab;
                   //TransformToLab(cellid0, sth->getPosition(), sLab);
                   sLab = sth->getPosition();
                   double pixelPos[3];
                   pixelPos[0] = sLab[0] / _pixelSizeX;
                   pixelPos[1] = sLab[1] / _pixelSizeY;
                   newsth->setPosition(pixelPos);
                   newsth->setEDep(charge); // in unit of electrons
                   newsth->setTime(sth->getTime());
                   newsth->setPathLength(simTrkHit->getPathLength());
                   newsth->setMCParticle(simTrkHit->getMCParticle());
                   newsth->setMomentum(simTrkHit->getMomentum());
                   newsth->setProducedBySecondary(simTrkHit->isProducedBySecondary());
                   newsth->setOverlay(simTrkHit->isOverlay());
                   STHLocCol->addElement(newsth);
                   recoHit->rawHits().push_back(newsth);
                }
              }
            }

            streamlog_out (DEBUG7) << "- number of pixels: " << recoHit->getRawHits().size() << std::endl;
            streamlog_out (DEBUG7) << "- MC particle p=" << std::sqrt(simTrkHit->getMomentum()[0]*simTrkHit->getMomentum()[0]+simTrkHit->getMomentum()[1]*simTrkHit->getMomentum()[1]+simTrkHit->getMomentum()[2]*simTrkHit->getMomentum()[2]) << std::endl;
            streamlog_out (DEBUG7) << "- isSecondary = " << simTrkHit->isProducedBySecondary() << ", isOverlay = " << simTrkHit->isOverlay() << std::endl;
            streamlog_out (DEBUG6) << "- List of constituents (pixels/strips):" << std::endl;
            for (size_t iH = 0; iH < recoHit->rawHits().size(); ++iH) {
                SimTrackerHit *hit = dynamic_cast<SimTrackerHit*>(recoHit->rawHits().at(iH));
                streamlog_out (DEBUG6) << "  - " << iH << ": Edep (e-) = " << hit->getEDep() << ", t (ns) =" << hit->getTime() << std::endl;
            }
            streamlog_out (DEBUG7) << "--------------------------------" << std::endl;

            THcol->addElement(recoHit);

            for (int k=0; k < int(simTrkHitVec.size()); ++k)
            {
                SimTrackerHit *hit = simTrkHitVec[k];
                delete hit;
            }

        }

        streamlog_out(DEBUG) << "Number of produced hits: " << THcol->getNumberOfElements()  << std::endl;

        //**************************************************************************
        // Add collection to event
        //**************************************************************************    

         evt->addCollection( THcol , _outputCollectionName.c_str() ) ;
         evt->addCollection( relCol , _colVTXRelation.c_str() ) ;
         if (_produceFullPattern != 0)
            evt->addCollection(STHLocCol, "VTXPixels");

    }

    streamlog_out(DEBUG9) << " Done processing event: " << evt->getEventNumber() 
        << "   in run:  " << evt->getRunNumber() << std::endl ;

    _nEvt ++ ;
}

void MuonCVXDDigitiser::check(LCEvent *evt)
{}

void MuonCVXDDigitiser::end()
{
    streamlog_out(DEBUG) << "   end called  " << std::endl;
    delete _fluctuate;
}

/** Function calculates local coordinates of the sim hit 
 * in the given ladder and local momentum of particle. 
 * Also returns module number and ladder number.
 * Local coordinate system within the ladder 
 * is defined as following :  <br> 
 *    - x axis lies in the ladder plane and orthogonal to the beam axis <br>
 *    - y axis lies in the ladder plane and parallel to the beam axis <br>
 *    - z axis is perpendicular to the ladder plane <br>
 * 
 */
void MuonCVXDDigitiser::FindLocalPosition(SimTrackerHit *hit, 
                                          double *localPosition,
                                          double *localDirection)
{
    // Use SurfaceManager to calculate local coordinates
    const int cellID0 = hit->getCellID0() ;
    SurfaceMap::const_iterator sI = _map->find( cellID0 ) ;
    const dd4hep::rec::ISurface* surf = sI->second ;
    Vector3D oldPos( hit->getPosition()[0], hit->getPosition()[1], hit->getPosition()[2] );
    // We need it?
    if ( ! surf->insideBounds( dd4hep::mm * oldPos ) ) {

        streamlog_out( DEBUG3 ) << "  hit at " << oldPos
                                << " is not on surface "
                                << *surf
                                << " distance: " << surf->distance(  dd4hep::mm * oldPos )
                                << std::endl;
      _currentLayer = -1;
      return;
    }    
    
    
    Vector2D lv = surf->globalToLocal( dd4hep::mm * oldPos  ) ;
    // Store local position in mm
    localPosition[0] = lv[0] / dd4hep::mm ;
    localPosition[1] = lv[1] / dd4hep::mm ;

    // Add also z ccordinate
    Vector3D origin( surf->origin()[0], surf->origin()[1], surf->origin()[2]);
    localPosition[2] = ( dd4hep::mm * oldPos - dd4hep::cm * origin ).dot( surf->normal() ) / dd4hep::mm;

    double Momentum[3];
    EVENT::MCParticle *mcp = hit->getMCParticle();
    for (int j = 0; j < 3; ++j) {
      if (mcp) {
        Momentum[j] = mcp->getMomentum()[j] * dd4hep::GeV;
      } else {
        Momentum[j] = hit->getMomentum()[j];
      }
    }

    // as default put electron's mass
    _currentParticleMass = 0.510e-3 * dd4hep::GeV;
    if (hit->getMCParticle())
        _currentParticleMass = std::max(hit->getMCParticle()->getMass() * dd4hep::GeV, _currentParticleMass);

    _currentParticleMomentum = sqrt(pow(Momentum[0], 2) + pow(Momentum[1], 2) 
                                    + pow(Momentum[2], 2));                   
                         
    localDirection[0] = Momentum * surf->u();
    localDirection[1] = Momentum * surf->v();
    localDirection[2] = Momentum * surf->normal();

    _currentPhi = _currentLadder * 2.0 * _layerHalfPhi[_currentLayer] + _layerPhiOffset[_currentLayer];

}

void MuonCVXDDigitiser::ProduceIonisationPoints(SimTrackerHit *hit)
{
    streamlog_out( DEBUG6 ) << "Creating Ionization Points" << std::endl;
    double pos[3] = {0,0,0};
    double dir[3] = {0,0,0};
    double entry[3];
    double exit[3];

    // hit and pos are in mm
    FindLocalPosition(hit, pos, dir);
    if ( _currentLayer == -1)
      return;
 
    entry[2] = -_layerHalfThickness[_currentLayer]; 
    exit[2] = _layerHalfThickness[_currentLayer];

    for (int i = 0; i < 2; ++i) {
        entry[i] = pos[i] + dir[i] * (entry[2] - pos[2]) / dir[2];
        exit[i]= pos[i] + dir[i] * (exit[2] - pos[2]) / dir[2];
    }

    for (int i = 0; i < 3; ++i) {
        _currentLocalPosition[i] = pos[i];
        _currentEntryPoint[i] = entry[i];
        _currentExitPoint[i] = exit[i];
    }

    streamlog_out( DEBUG5 ) << "local position: " << _currentLocalPosition[0] << ", " << _currentLocalPosition[1] << ", " << _currentLocalPosition[2] << std::endl;

    double tanx = dir[0] / dir[2];
    double tany = dir[1] / dir[2];  
    
    // trackLength is in mm -> limit length at 1cm
    double trackLength = std::min(_maxTrkLen,
         _layerThickness[_currentLayer] * sqrt(1.0 + pow(tanx, 2) + pow(tany, 2)));
  
    _numberOfSegments = ceil(trackLength / _segmentLength );
    double dEmean = (dd4hep::keV * _energyLoss * trackLength) / ((double)_numberOfSegments);
    _ionisationPoints.resize(_numberOfSegments);

    streamlog_out( DEBUG6 ) <<  "Track path length: " << trackLength << ", calculated dEmean * N_segment = " << dEmean << " * " << _numberOfSegments << " = " << dEmean*_numberOfSegments << std::endl;
    _eSum = 0.0;

    // TODO _segmentLength may be different from segmentLength, is it ok?
    double segmentLength = trackLength / ((double)_numberOfSegments);
    _segmentDepth = _layerThickness[_currentLayer] / ((double)_numberOfSegments);

    double z = -_layerHalfThickness[_currentLayer] - 0.5 * _segmentDepth;
    
    double hcharge = ( hit->getEDep() / dd4hep::GeV ); 	
    streamlog_out (DEBUG5) << "Number of ionization points: " << _numberOfSegments << ", G4 EDep = "  << hcharge << std::endl;
    for (int i = 0; i < _numberOfSegments; ++i)
    {
        z += _segmentDepth;
        double x = pos[0] + tanx * (z - pos[2]);
        double y = pos[1] + tany * (z - pos[2]);
        // momentum in MeV/c, mass in MeV, tmax (delta cut) in MeV, 
        // length in mm, meanLoss eloss in MeV.
        double de = _fluctuate->SampleFluctuations(double(_currentParticleMomentum * dd4hep::keV / dd4hep::MeV),
                                                   double(_currentParticleMass * dd4hep::keV / dd4hep::MeV),
                                                   _cutOnDeltaRays,
                                                   segmentLength,
                                                   double(dEmean / dd4hep::MeV)) * dd4hep::MeV;
        _eSum = _eSum + de;

        IonisationPoint ipoint;
        ipoint.eloss = de;
        ipoint.x = x;
        ipoint.y = y;
        ipoint.z = z;
        _ionisationPoints[i] = ipoint;
        streamlog_out (DEBUG2) << " " << i << ": z=" << z << ", eloss = " << de << "(total so far: " << _eSum << "), x=" << x << ", y=" << y << std::endl;
    }
   
    const double thr = _deltaEne/_electronsPerKeV * dd4hep::keV;
    while ( hcharge > _eSum + thr ) {
      // Add additional charge sampled from an 1 / n^2 distribution.
      // Adjust charge to match expectations
      const double       q = randomTail( thr, hcharge - _eSum );
      const unsigned int h = floor(RandFlat::shoot(0.0, (double)_numberOfSegments ));
      _ionisationPoints[h].eloss += q;
      _eSum += q;
    }
    streamlog_out (DEBUG5) << "Padding each segment charge (1/n^2 pdf) until total below " << _deltaEne << "e- threshold. New total energy: " << _eSum << std::endl;
    streamlog_out (DEBUG3) << "List of ionization points:" << std::endl;
    for (int i =0; i < _numberOfSegments; ++i) {
        streamlog_out (DEBUG3) << "- " << i << ": E=" << _ionisationPoints[i].eloss 
            << ", x=" << _ionisationPoints[i].x << ", y=" << _ionisationPoints[i].y << ", z=" << _ionisationPoints[i].z << std::endl;
    }
}

void MuonCVXDDigitiser::ProduceSignalPoints()
{
    _signalPoints.resize(_numberOfSegments);

    // run over ionisation points
    streamlog_out (DEBUG6) << "Creating signal points" << std::endl;
    for (int i = 0; i < _numberOfSegments; ++i)
    {
        IonisationPoint ipoint = _ionisationPoints[i];
        double z = ipoint.z;
        double x = ipoint.x;
        double y = ipoint.y;
        double DistanceToPlane = _layerHalfThickness[_currentLayer] - z;
        double xOnPlane = x + _tanLorentzAngleX * DistanceToPlane;
        double yOnPlane = y + _tanLorentzAngleY * DistanceToPlane;

        // For diffusion-coeffieint calculation, see e.g. https://www.slac.stanford.edu/econf/C060717/papers/L008.PDF
        // or directly Eq. 13 of https://cds.cern.ch/record/2161627/files/ieee-tns-07272141.pdf
        // diffusionCoefficient = sqrt(2*D / mu / V), where
        //  - D = 12 cm^2/s // diffusion constant
        //  - mu = 450 cm^2/s/V // mobility
        //  - V = 10-30 V // expected depletion voltage
        //  => _diffusionCoefficient = 0.04-0.07
        // and diffusion sigma = _diffusionCoefficient * DistanceToPlane
        // e.g. fot 50um diffusion sigma = 2.1 - 3.7 um
        //double DriftLength = DistanceToPlane * sqrt(1.0 + pow(_tanLorentzAngleX, 2) 
        //                                                + pow(_tanLorentzAngleY, 2));
        double SigmaDiff = DistanceToPlane * _diffusionCoefficient;

        double SigmaX = SigmaDiff * sqrt(1.0 + pow(_tanLorentzAngleX, 2));
        double SigmaY = SigmaDiff * sqrt(1.0 + pow(_tanLorentzAngleY, 2));

        // energy is in keV       
        double charge = (ipoint.eloss / dd4hep::keV) * _electronsPerKeV;

        SignalPoint  spoint;
        spoint.x = xOnPlane;
        spoint.y = yOnPlane;
        spoint.sigmaX = SigmaX;
        spoint.sigmaY = SigmaY;
        spoint.charge = charge; // electrons x keV
        _signalPoints[i] = spoint;
	    streamlog_out (DEBUG3) << "- " << i << ": charge=" << charge 
            << ", x="<<xOnPlane << "(delta=" << xOnPlane - x << ")"
            << ", y="<<yOnPlane << "(delta=" << yOnPlane - y << ")"
            << ", sigmaDiff=" << SigmaDiff
            << ", sigmaX="<<SigmaX <<", sigmay="<<SigmaY << std::endl;
    }
}

void MuonCVXDDigitiser::ProduceHits(SimTrackerHitImplVec &simTrkVec, SimTrackerHit &simHit)
{  
    simTrkVec.clear();
    std::map<int, SimTrackerHitImpl*> hit_Dict;
    streamlog_out (DEBUG6) << "Creating hits" << std::endl;

    for (int i=0; i<_numberOfSegments; ++i)
    {
        SignalPoint spoint = _signalPoints[i];
        double xCentre = spoint.x;
        double yCentre = spoint.y;
        double sigmaX = spoint.sigmaX;
        double sigmaY = spoint.sigmaY;
        double xLo = spoint.x - 3 * spoint.sigmaX;
        double xUp = spoint.x + 3 * spoint.sigmaX;
        double yLo = spoint.y - 3 * spoint.sigmaY;
        double yUp = spoint.y + 3 * spoint.sigmaY;
        
        int ixLo, ixUp, iyLo, iyUp;

        TransformXYToCellID(xLo, yLo, ixLo, iyLo);
        TransformXYToCellID(xUp, yUp, ixUp, iyUp);
        streamlog_out (DEBUG5) << i << ": Pixel idx boundaries: ixLo=" << ixLo << ", iyLo=" << iyLo  
            <<  ", ixUp=" << ixUp << ", iyUp=" << iyUp << std::endl;

        for (int ix = ixLo; ix< ixUp + 1; ++ix)
        {
            if ( (ix < 0) or (ix >= GetPixelsInaColumn()) ) continue;

            for (int iy = iyLo; iy < iyUp + 1; ++iy)
            {
                if ( (iy < 0) or (iy >= GetPixelsInaRow()) ) continue;

                double xCurrent, yCurrent;
                TransformCellIDToXY(ix, iy, xCurrent, yCurrent);
                
                gsl_sf_result result;
                /*int status = */gsl_sf_erf_Q_e((xCurrent - 0.5 * _pixelSizeX - xCentre)/sigmaX, &result);
                double LowerBound = 1 - result.val;

                /*status = */gsl_sf_erf_Q_e((xCurrent + 0.5 * _pixelSizeX - xCentre)/sigmaX, &result);
                double UpperBound = 1 - result.val;
                double integralX = UpperBound - LowerBound;

                /*status = */gsl_sf_erf_Q_e((yCurrent - 0.5 * _pixelSizeY - yCentre)/sigmaY, &result);
                LowerBound = 1 - result.val;

                /*status = */gsl_sf_erf_Q_e((yCurrent + 0.5 * _pixelSizeY - yCentre)/sigmaY, &result);
                UpperBound = 1 - result.val;
                double integralY = UpperBound - LowerBound;

                float totCharge = float(spoint.charge * integralX * integralY);

                int pixelID = GetPixelsInaRow() * ix + iy;
              
                auto item = hit_Dict.find(pixelID);
                if (item == hit_Dict.end())
                {
                    SimTrackerHitImpl *tmp_hit = new SimTrackerHitImpl();
                    double pos[3] = {
                        xCurrent,
                        yCurrent,
                        _layerHalfThickness[_currentLayer]
                    };
                    tmp_hit->setPosition(pos);
                    tmp_hit->setCellID0(pixelID);                   // workaround: cellID used for pixel index
                    tmp_hit->setEDep(totCharge);
                    tmp_hit->setTime(simHit.getTime()); //usual true timing as starting point
                  
                    hit_Dict.emplace(pixelID, tmp_hit);
                    streamlog_out (DEBUG2) << "Created new pixel hit at idx=" << ix << ", idy=" << iy << ", charge=" << totCharge << std::endl;
                }
                else
                {
                    float edep = item->second->getEDep();
                    edep += totCharge;
                    item->second->setEDep(edep);
                    //TODO: handle multiple times. For now not needed since all deposits arrive at the same true time.
                    streamlog_out (DEBUG2) << "Updating pixel hit at idx=" << ix << ", idy=" << iy << ", total charge=" << edep << "(delta = " << totCharge << ")" << std::endl;
                }
            }
        }
    }
    streamlog_out (DEBUG4) << "List of pixel hits created:" << std::endl;
    int idx=0;
    for(auto item : hit_Dict)
    {
        simTrkVec.push_back(item.second);       
        streamlog_out (DEBUG4) << idx++ << ": x=" << item.second->getPosition()[0] << ", y=" << item.second->getPosition()[1] << ", z=" << item.second->getPosition()[2] << ", EDep = " << item.second->getEDep() << std::endl;
    }
}

/**
 * Function that fluctuates charge (in units of electrons)
 * deposited on the fired pixels according to the Poisson
 * distribution...
 */
void MuonCVXDDigitiser::PoissonSmearer(SimTrackerHitImplVec &simTrkVec)
{
    streamlog_out (DEBUG6) << "Adding Poisson smear to charge" << std::endl;
    for (int ihit = 0; ihit < int(simTrkVec.size()); ++ihit)
    {
        SimTrackerHitImpl *hit = simTrkVec[ihit];
        float charge = hit->getEDep();
        float rng;
        if (charge > 1e+03) // assume Gaussian
        {
            rng = float(RandGauss::shoot(charge, sqrt(charge)));
        }
        else // assume Poisson
        {
            rng = float(RandPoisson::shoot(charge));
        }
        hit->setEDep(rng);
        streamlog_out (DEBUG4) << ihit << ": x=" << hit->getPosition()[0] << ", y=" << hit->getPosition()[1] << ", z=" << hit->getPosition()[2] 
            << ", charge = " << rng << "(delta = " << charge-rng << ")" << std::endl;
    }  
}

/**
 * Simulation of electronic noise.
 */
void MuonCVXDDigitiser::GainSmearer(SimTrackerHitImplVec &simTrkVec)
{
    streamlog_out (DEBUG6) << "Adding FE noise smear to charge" << std::endl;
    for (int i = 0; i < (int)simTrkVec.size(); ++i)
    {
        double Noise = RandGauss::shoot(0., _electronicNoise);
        SimTrackerHitImpl *hit = simTrkVec[i];
        hit->setEDep(hit->getEDep() + float(Noise));
        streamlog_out (DEBUG4) << i << ": x=" << hit->getPosition()[0] << ", y=" << hit->getPosition()[1] << ", z=" << hit->getPosition()[2] 
            << ", charge = " << hit->getEDep() << "(delta = " << Noise << ")" << std::endl;
    }
}

/**
 * Apply threshold.  
 * Sets the charge to 0 if less than the threshold
 * Smears the threshold by a Gaussian if sigma > 0
 */
void MuonCVXDDigitiser::ApplyThreshold(SimTrackerHitImplVec &simTrkVec)
{
   streamlog_out (DEBUG6) << "Applying threshold" << std::endl;
   float actualThreshold = _threshold;
   
   for (int i = 0; i < (int)simTrkVec.size(); ++i)
   {
     SimTrackerHitImpl *hit = simTrkVec[i];
     
     double smear = 0;
     float origCharge = hit->getEDep();

     if (_thresholdSmearSigma > 0) smear = RandGauss::shoot(0., _thresholdSmearSigma);

     actualThreshold = actualThreshold + smear;
     if (hit->getEDep() <= actualThreshold) hit->setEDep(0.0);
     
     streamlog_out (DEBUG4) << i << ": x=" << hit->getPosition()[0] << ", y=" << hit->getPosition()[1] << ", z=" << hit->getPosition()[2]
			   << ", new charge = " << hit->getEDep() << ", previous charge = " << origCharge
			   << " smeared threshold = " << actualThreshold << "(delta = " << smear << ")" << std::endl;
   }
}

/**
 * Digitizes the charge.
 * Discretization based on number of bits and bin width scheme.
 */
void MuonCVXDDigitiser::ChargeDigitizer(SimTrackerHitImplVec &simTrkVec)
{
  streamlog_out (DEBUG6) << "Charge discretization" << std::endl;
  
  float minThreshold = _threshold;
  float maxThreshold = _chargeMax;
  //int split = 0.3; -- future use
  int numBins = pow(2, _ChargeDigitizeNumBits)-1;

  double discCharge=-999; 

  for (int i = 0; i < (int)simTrkVec.size(); ++i) {
    SimTrackerHitImpl *hit = simTrkVec[i];

    float origCharge =	hit->getEDep();
    discCharge = origCharge;
     
    switch(_ChargeDigitizeBinning) {
        case 0: { // uniform binning
            if (origCharge < 1.0) break;
    	    float binWidth = (maxThreshold-minThreshold)/(numBins);
	        if (origCharge < binWidth) discCharge = (minThreshold+binWidth)/2;
	        else if (origCharge > maxThreshold) discCharge = (maxThreshold-binWidth/2);
            else discCharge = ((ceil((origCharge-binWidth)/binWidth)*binWidth)*2+binWidth)/2;
	        break;
        }
        case 1: { // variable binning
	        if (origCharge < 1.0) break;
	        int binVal=-1;
	        for(unsigned int idx = 0; idx < _DigitizedBins.size()-1; idx++) {
	            if (_DigitizedBins[idx+1] > origCharge) {
		            binVal = idx;
		            break;
	            }
            }
	        if (binVal < 0) discCharge = (_DigitizedBins[_DigitizedBins.size()-2] + _DigitizedBins[_DigitizedBins.size()-1]) / 2;
	        else discCharge = (_DigitizedBins[binVal] + _DigitizedBins[binVal+1]) / 2;
	        break;
	    }
    }
    hit->setEDep(discCharge);

    streamlog_out (DEBUG4) << i << ": x=" << hit->getPosition()[0] << ", y=" << hit->getPosition()[1] << ", z=" << hit->getPosition()[2]
                        << ", new charge = " << hit->getEDep() << ", previous charge = " << origCharge
			            << ", number of bits = " << _ChargeDigitizeNumBits
			            << ", binning scheme = " << _ChargeDigitizeBinning << std::endl;
    }
}

/**
 * Apply effective measurement resolution.
 * TODO: Right now assuming completely uncorrelated resolution across pixels, will need to divide into:
 * - correlated across pixels, uncorrelated across clusters
 * - correlated within the event, un-correlate 
*/
void MuonCVXDDigitiser::TimeSmearer(SimTrackerHitImplVec &simTrkVec)
{
    streamlog_out (DEBUG6) << "Adding resolution effect to timing measurements" << std::endl;
    for (int i = 0; i < (int)simTrkVec.size(); ++i)
    {
        float delta = RandGauss::shoot(0., _timeSmearingSigma);
        SimTrackerHitImpl *hit = simTrkVec[i];
        hit->setTime(hit->getTime() + delta);
        streamlog_out (DEBUG4) << i << ": x=" << hit->getPosition()[0] << ", y=" << hit->getPosition()[1] << ", z=" << hit->getPosition()[2] 
            << ", time = " << hit->getTime() << "(delta = " << delta << ")" << std::endl;
    }
}

/**
 * Digitizes the time information.
 * Discretization based on number of bits and bin width scheme.
 */
void MuonCVXDDigitiser::TimeDigitizer(SimTrackerHitImplVec &simTrkVec)
{
    streamlog_out (DEBUG6) << "Time discretization" << std::endl;
  
    static const int numBins = pow(2, _TimeDigitizeNumBits)-1;

    double discTime;

    for (int i = 0; i < (int)simTrkVec.size(); ++i)
    {
        SimTrackerHitImpl *hit = simTrkVec[i];

        float origTime = hit->getTime();
        discTime = origTime;
     
        switch(_TimeDigitizeBinning) {
        case 0: // uniform binning
	        static const float binWidth = _timeMax/numBins;
	        if (origTime < binWidth) discTime = binWidth/2;
	        else if (origTime > _timeMax) discTime = _timeMax - binWidth/2;
            else discTime = ((ceil((origTime-binWidth)/binWidth)*binWidth)*2+binWidth)/2;
	        break;
        default:
            streamlog_out(ERROR) << "Invalid setting for pixel time digitization binning. Retaining original time." << std::endl;
        }
        hit->setTime(discTime);

        streamlog_out (DEBUG4) << i << ": x=" << hit->getPosition()[0] << ", y=" << hit->getPosition()[1] << ", z=" << hit->getPosition()[2]
                << ", new time = " << hit->getTime() << ", previous time = " << origTime << std::endl;
    } //end loop over pixel cells
}

/**
 * Emulates reconstruction of Tracker Hit 
 * Tracker hit position is reconstructed as weighted average of edge pixels.
 * The position is corrected for Lorentz shift.
 * Time is the arithmetic average of constituents.
 */
TrackerHitPlaneImpl *MuonCVXDDigitiser::ReconstructTrackerHit(SimTrackerHitImplVec &simTrkVec)
{
    double pos[3] = {0, 0, 0};

    double minX = 99999999;
    double maxX = -99999999;
    std::vector<int> minIdx;
    std::vector<int> maxIdx;
      
    double charge = 0; // total cluster charge
    unsigned int size = 0; // number of pixels in the cluster
    unsigned int edge_size = 0; //number of pixels at the edge of cluster
    streamlog_out (DEBUG6) << "Creating reconstructed cluster" << std::endl;
    double time = 0; //average time

    /* Get extreme positions, currently only implemented for barrel */
    /* Calculate the mean */
    for (size_t iHit=0; iHit < simTrkVec.size(); ++iHit)
    {
        SimTrackerHit *hit = simTrkVec[iHit];

        //check for non-zero value (pixels below threshold have already been set to zero)
	    if (hit->getEDep() < 1.0) continue;
	
        size += 1;
        time += hit->getTime();
        charge += hit->getEDep();
        pos[0] += hit->getPosition()[0];
        pos[1] += hit->getPosition()[1];
        streamlog_out (DEBUG0) << iHit << ": Averaging position, x=" << hit->getPosition()[0] << ", y=" << hit->getPosition()[1] << ", weight(EDep)=" << hit->getEDep() << std::endl;

        if (hit->getPosition()[0] < minX) {
	        minX = hit->getPosition()[0];
	        minIdx.push_back(iHit);
	    }
	    if (hit->getPosition()[0] > maxX) {
          maxX = hit->getPosition()[0];
          maxIdx.push_back(iHit);
	    }
    }
    
    for (unsigned int i = 0; i < minIdx.size(); i++){
      SimTrackerHit *hit = simTrkVec[minIdx[i]];
      
      if (minX < hit->getPosition()[0]) continue;
      streamlog_out (DEBUG4) << "min val = " << minX << ", hit position = " <<  hit->getPosition()[0]  << std::endl;
      edge_size += 1;
      pos[0] += hit->getPosition()[0];
      pos[1] += hit->getPosition()[1];
    }

    for (unsigned int i = 0; i < maxIdx.size(); i++){
      SimTrackerHit *hit = simTrkVec[maxIdx[i]];

      if (maxX > hit->getPosition()[0]) continue;
      streamlog_out (DEBUG4) << "max val = " << maxX << ", hit position = " <<  hit->getPosition()[0]  << std::endl;
      edge_size += 1;
      pos[0] += hit->getPosition()[0];
      pos[1] += hit->getPosition()[1];
    }

    if ( not (charge > 0.) ) return nullptr;

    TrackerHitPlaneImpl *recoHit = new TrackerHitPlaneImpl();
    recoHit->setEDep((charge / _electronsPerKeV) * dd4hep::keV);

    pos[0] /= edge_size;
    streamlog_out (DEBUG1) << "Position: x = " << pos[0] << " + " << _layerHalfThickness[_currentLayer] * _tanLorentzAngleX << "(LA-correction)";
    pos[0] -= _layerHalfThickness[_currentLayer] * _tanLorentzAngleX;
    streamlog_out (DEBUG1) << " = " << pos[0];
    pos[1] /= edge_size;
    streamlog_out (DEBUG1) << "; y = " << pos[1] << " + " << _layerHalfThickness[_currentLayer] * _tanLorentzAngleY << "(LA-correction)";
    pos[1] -= _layerHalfThickness[_currentLayer] * _tanLorentzAngleY;
    streamlog_out (DEBUG1) << " = " << pos[1];

    recoHit->setPosition(pos);
    recoHit->setdU( _pixelSizeX / sqrt(12) );
    recoHit->setdV( _pixelSizeY / sqrt(12) );

    time /= size;
    recoHit->setTime(time);
    streamlog_out (DEBUG1) << ", time (ns) = " << time << std::endl;
          
    return recoHit;
}

/** Function transforms local coordinates in the ladder
 * into global coordinates
 */
void MuonCVXDDigitiser::TransformToLab(const int cellID, const double *xLoc, double *xLab)
{
    // Use SurfaceManager to calculate global coordinates
    SurfaceMap::const_iterator sI = _map->find( cellID ) ;
    const dd4hep::rec::ISurface* surf = sI->second ;
    Vector2D oldPos( xLoc[0] * dd4hep::mm, xLoc[1] * dd4hep::mm );
    Vector3D lv = surf->localToGlobal( oldPos  ) ;
    // Store local position in mm
    for ( int i = 0; i < 3; i++ )
      xLab[i] = lv[i] / dd4hep::mm;
}

/**
 * Function calculates position in pixel matrix based on the 
 * local coordinates of point in the ladder.
 */
void MuonCVXDDigitiser::TransformXYToCellID(double x, double y, int & ix, int & iy)
{
    int layer = _currentLayer;

    // Shift all of L/2 so that all numbers are positive
    double yInLadder = y + _layerLadderLength[layer] / 2;
    iy = int(yInLadder / _pixelSizeY);

    double xInLadder = x + _layerLadderHalfWidth[layer];
    ix = int(xInLadder / _pixelSizeX);
}

/**
 Function calculates position in the local frame 
 based on the index of pixel in the ladder.
*/
void MuonCVXDDigitiser::TransformCellIDToXY(int ix, int iy, double & x, double & y)
{
    int layer = _currentLayer;
    // Put the point in the cell center
    y = ((0.5 + double(iy)) * _pixelSizeY) - _layerLadderLength[layer] / 2;
    x = ((0.5 + double(ix)) * _pixelSizeX) - _layerLadderHalfWidth[layer];
}

int MuonCVXDDigitiser::GetPixelsInaColumn() //SP: why columns!?! I would have guess row..
{
    return ceil(_layerLadderWidth[_currentLayer] / _pixelSizeX);
}

int MuonCVXDDigitiser::GetPixelsInaRow()
{
    return ceil(_layerLadderLength[_currentLayer] / _pixelSizeY);
}

void MuonCVXDDigitiser::PrintGeometryInfo()
{
    streamlog_out(MESSAGE) << "Number of layers: " << _numberOfLayers << std::endl;
    streamlog_out(MESSAGE) << "Pixel size X: " << _pixelSizeX << std::endl;
    streamlog_out(MESSAGE) << "Pixel size Y: " << _pixelSizeY << std::endl;
    streamlog_out(MESSAGE) << "Electrons per KeV: " << _electronsPerKeV << std::endl;
    streamlog_out(MESSAGE) << "Segment depth: " << _segmentDepth << std::endl;
    for (int i = 0; i < _numberOfLayers; ++i) 
    {
        streamlog_out(MESSAGE) << "Layer " << i << std::endl;
        streamlog_out(MESSAGE) << "  Number of ladders: " << _laddersInLayer[i] << std::endl;
        streamlog_out(MESSAGE) << "  Radius: " << _layerRadius[i] << std::endl;
        streamlog_out(MESSAGE) << "  Ladder length: " << _layerLadderLength[i] << std::endl;
        streamlog_out(MESSAGE) << "  Ladder width: "<< _layerLadderWidth[i] << std::endl;
        streamlog_out(MESSAGE) << "  Ladder half width: " << _layerLadderHalfWidth[i] << std::endl;
        streamlog_out(MESSAGE) << "  Phi offset: " << _layerPhiOffset[i] << std::endl;
        streamlog_out(MESSAGE) << "  Active Si offset: " << _layerActiveSiOffset[i] << std::endl;
        streamlog_out(MESSAGE) << "  Half phi: " << _layerHalfPhi[i] << std::endl;
        streamlog_out(MESSAGE) << "  Thickness: " << _layerThickness[i] << std::endl;
        streamlog_out(MESSAGE) << "  Half thickness: " << _layerHalfThickness[i] << std::endl;
    }
}


//=============================================================================
// Sample charge from 1 / n^2 distribution.
//=============================================================================
double MuonCVXDDigitiser::randomTail( const double qmin, const double qmax ) {

  const double offset = 1. / qmax;
  const double range  = ( 1. / qmin ) - offset;
  const double u      = offset + RandFlat::shoot() * range;
  return 1. / u;
}



