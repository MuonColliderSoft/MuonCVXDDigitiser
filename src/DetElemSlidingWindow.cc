#include "DetElemSlidingWindow.h"
#include <EVENT/MCParticle.h>
#include "DDRec/DetectorData.h"
#include "DD4hep/DD4hepUnits.h"
#include "marlin/VerbosityLevels.h"
#include "streamlog/streamlog.h"

#include "gsl/gsl_sf_erf.h"
#include "gsl/gsl_math.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandPoisson.h"
#include "CLHEP/Random/RandFlat.h"

#include <iostream>

using std::max;
using std::min;
using std::sqrt;
using std::atan;
using std::vector;
using dd4hep::rec::ISurface;
using dd4hep::rec::Vector2D;
using dd4hep::rec::Vector3D;
using CLHEP::RandGauss;
using CLHEP::RandPoisson;
using CLHEP::RandFlat;

DetElemSlidingWindow::DetElemSlidingWindow(HitTemporalIndexes& htable,
                                           AbstractSensor& sensor,
                                           float wsize,
                                           float starttime,
                                           double tanLorentzAngleX,
                                           double tanLorentzAngleY,
                                           double cutOnDeltaRays,
                                           double diffusionCoefficient,
                                           double electronsPerKeV,
                                           double segmentLength,
                                           double energyLoss,
                                           double widthOfCluster,
                                           double electronicNoise,
                                           double maxTrkLen,
                                           double maxEnergyDelta,
                                           const SurfaceMap* s_map,
                                           bool zSegmented,
                                           bool isBarrel):
    curr_time(starttime + wsize / 2),  // window centered in the middle
    time_click(wsize),
    _htable(htable),
    _sensor(sensor),
    _tanLorentzAngleX(tanLorentzAngleX),
    _tanLorentzAngleY(tanLorentzAngleY),
    _cutOnDeltaRays(cutOnDeltaRays),
    _diffusionCoefficient(diffusionCoefficient),
    _electronsPerKeV(electronsPerKeV),
    _segmentLength(segmentLength),
    _energyLoss(energyLoss),
    _widthOfCluster(widthOfCluster),
    _electronicNoise(electronicNoise),
    _maxTrkLen(maxTrkLen),
    _deltaEne(maxEnergyDelta),
    signals(),
    surf_map(s_map),
    _zSegmented(zSegmented),
    _isBarrel(isBarrel),
    cell_decoder(sensor.GetCellIDFormatStr())
{
    _fluctuate = new G4UniversalFluctuation();
}

DetElemSlidingWindow::~DetElemSlidingWindow()
{
    delete(_fluctuate);
}


bool DetElemSlidingWindow::active()
{
    bool hasMoreHits = _htable.GetHitNumber(_sensor.GetLayer(), _sensor.GetLadder()) > 0;
    bool sensorOn = _sensor.IsActive();
    return hasMoreHits || sensorOn;
}

int DetElemSlidingWindow::process()
{
    float window_radius = time_click / 2;

    for (SimTrackerHit* hit = _htable.CurrentHit(_sensor.GetLayer(), _sensor.GetLadder());
         hit != nullptr && hit->getTime() - curr_time < window_radius;
         hit = _htable.CurrentHit(_sensor.GetLayer(), _sensor.GetLadder()))
    {
        if (streamlog::out.write<streamlog::DEBUG6>())
#pragma omp critical
        {
            float mcp_r = sqrt(pow(hit->getPosition()[0], 2) + pow(hit->getPosition()[1], 2));
            float mcp_phi = atan(hit->getPosition()[1] / hit->getPosition()[0]);
            float mcp_theta = hit->getPosition()[2] == 0 ? 3.1416/2 : atan(mcp_r / hit->getPosition()[2]);
            double mom_norm = sqrt(pow(hit->getMomentum()[0], 2) + pow(hit->getMomentum()[1], 2)
                                   + pow(hit->getMomentum()[2], 2));
            int segment_id = cell_decoder(hit)["sensor"];
            streamlog::out() << "Processing simHit from layer = " << _sensor.GetLayer()
                             << ", ladder = " << _sensor.GetLadder() 
                             << ", sensor = " << segment_id << std::endl
                             << "Time window centered in " << curr_time
                             << ", Hits available = "
                             << _htable.GetHitNumber(_sensor.GetLayer(), _sensor.GetLadder())
                             << std::endl
                             << "- EDep = " << hit->getEDep() * dd4hep::GeV / dd4hep::keV
                             << " keV, path length = " << hit->getPathLength() * 1000.
                             << " um" << std::endl
                             << "- Position (mm) x,y,z,t = " << hit->getPosition()[0] << ", "
                             << hit->getPosition()[1] << ", " << hit->getPosition()[2]
                             << ", " << hit->getTime() << std::endl
                             << "- Position r(mm),phi,theta = " << mcp_r << ", " << mcp_phi
                             << ", " << mcp_theta << std::endl
                             << "- MC particle pdg = " << hit->getMCParticle()->getPDG() << std::endl
                             << "- MC particle p (GeV) = " << mom_norm << std::endl
                             << "- isSecondary = " << hit->isProducedBySecondary()
                             << ", isOverlay = " << hit->isOverlay() << std::endl
                             << "- Quality = " << hit->getQuality() << std::endl;
        }

        StoreSignalPoints(hit);
        _htable.DisposeHit(_sensor.GetLayer(), _sensor.GetLadder());
    }

    if (!signals.empty())
    {
        streamlog_out(DEBUG) << "Signal points for " << _sensor.GetLayer() << ":" << _sensor.GetLadder()
                               << " = " << signals.size() << std::endl;

        for (TimedSignalPoint spoint = signals.front();
             curr_time - spoint.sim_hit->getTime() > window_radius;
             spoint = signals.front())
        {
            signals.pop_front();
            if (signals.empty()) break;
        }
    }

    UpdatePixels();
    curr_time += time_click;

    return signals.size();
}

float DetElemSlidingWindow::get_time()
{
    return curr_time;
}

void DetElemSlidingWindow::UpdatePixels()
{
    _sensor.InitHitRegister();
    _sensor.BeginClockStep();

    float window_radius = time_click / 2;

    for (auto spoint : signals)
    {
        if (spoint.sim_hit->getTime() > curr_time + window_radius) break;

        double xHFrame = _widthOfCluster * spoint.sigmaX;
        double yHFrame = _widthOfCluster * spoint.sigmaY;
        
        int ixLo = _sensor.XToPixelRow(spoint.x - xHFrame);
        int iyLo = _sensor.YToPixelCol(spoint.y - yHFrame);

        int ixUp = _sensor.XToPixelRow(spoint.x + xHFrame);
        int iyUp = _sensor.YToPixelCol(spoint.y + yHFrame);

        for (int ix = ixLo; ix < ixUp + 1; ++ix)
        {
            if (ix < 0 || ix >= _sensor.GetLadderRows()) continue;

            for (int iy = iyLo; iy < iyUp + 1; ++iy)
            {
                if (iy < 0 || iy >= _sensor.GetLadderCols()) continue;

                double xCurrent = _sensor.PixelRowToX(ix);
                double yCurrent = _sensor.PixelColToY(iy);
                
                gsl_sf_result result;
                int status = gsl_sf_erf_Q_e((xCurrent - 0.5 * _sensor.GetPixelSizeX() - spoint.x) / spoint.sigmaX, &result);
                double LowerBound = 1 - result.val;

                status = gsl_sf_erf_Q_e((xCurrent + 0.5 * _sensor.GetPixelSizeX() - spoint.x) / spoint.sigmaX, &result);
                double UpperBound = 1 - result.val;
                double integralX = UpperBound - LowerBound;

                status = gsl_sf_erf_Q_e((yCurrent - 0.5 * _sensor.GetPixelSizeY() - spoint.y) / spoint.sigmaY, &result);
                LowerBound = 1 - result.val;

                status = gsl_sf_erf_Q_e((yCurrent + 0.5 * _sensor.GetPixelSizeY() - spoint.y) / spoint.sigmaY, &result);
                UpperBound = 1 - result.val;
                double integralY = UpperBound - LowerBound;

                float totCharge = float(spoint.charge * integralX * integralY);

                _sensor.UpdatePixel(ix, iy, totCharge);
                _sensor.RegisterHit(ix, iy, spoint.sim_hit);
            }
        }
    }

    _sensor.EndClockStep();
}

void DetElemSlidingWindow::StoreSignalPoints(SimTrackerHit* hit)
{
    // hit and pos are in mm
    double pos[3] = {0,0,0};
    double dir[3] = {0,0,0};
    double entry[3];
    double exit[3];

    // ************************* Find local position **************************
    SurfaceMap::const_iterator sI = surf_map->find(hit->getCellID0()) ;
    const ISurface* surf = sI->second ;

    Vector3D oldPos( hit->getPosition()[0], hit->getPosition()[1], hit->getPosition()[2] );

    if (!surf->insideBounds(dd4hep::mm * oldPos))
    {
        if (streamlog::out.write<streamlog::DEBUG6>())
#pragma omp critical
        {
            streamlog::out() << "  hit at " << oldPos << " is not on surface " << *surf
                             << " distance: " << surf->distance(dd4hep::mm * oldPos) << std::endl;
        }
        return;
    }    


    Vector2D lv = surf->globalToLocal( dd4hep::mm * oldPos  ) ;
    // Store local position in mm
    pos[0] = lv[0] / dd4hep::mm;
    pos[1] = lv[1] / dd4hep::mm;
    if (_isBarrel && _zSegmented) {
        // See MuonCVXDDigitiser::processEvent
        int segment_id = cell_decoder(hit)["sensor"];

        float s_offset = _sensor.GetSensorCols() * _sensor.GetPixelSizeY() * (float(segment_id) + 0.5);
        s_offset -= _sensor.GetHalfLength();
        pos[1] += s_offset;
    }

    // Add also z ccordinate
    Vector3D origin( surf->origin()[0], surf->origin()[1], surf->origin()[2]);
    pos[2] = ( dd4hep::mm * oldPos - dd4hep::cm * origin ).dot( surf->normal() ) / dd4hep::mm;

    double Momentum[3];
    for (int j = 0; j < 3; ++j) 
      if (hit->getMCParticle())
        Momentum[j] = hit->getMCParticle()->getMomentum()[j] * dd4hep::GeV;
      else
        Momentum[j] = hit->getMomentum()[j];

    // as default put electron's mass
    double particleMass = 0.510e-3 * dd4hep::GeV;
    if (hit->getMCParticle())
        particleMass = max(hit->getMCParticle()->getMass() * dd4hep::GeV, particleMass);

    double particleMomentum = sqrt(pow(Momentum[0], 2) + pow(Momentum[1], 2) + pow(Momentum[2], 2));                   
                         
    dir[0] = Momentum * surf->u();
    dir[1] = Momentum * surf->v();
    dir[2] = Momentum * surf->normal();

    // ************************************************************************

    entry[2] = -_sensor.GetHalfThickness(); 
    exit[2] = _sensor.GetHalfThickness();

    for (int i = 0; i < 2; ++i) {
        entry[i] = pos[i] + dir[i] * (entry[2] - pos[2]) / dir[2];
        exit[i]= pos[i] + dir[i] * (exit[2] - pos[2]) / dir[2];
    }

    double tanx = dir[0] / dir[2];
    double tany = dir[1] / dir[2];  
    
    // trackLength is in mm
    double trackLength = min(_maxTrkLen, _sensor.GetThickness() * sqrt(1.0 + pow(tanx, 2) + pow(tany, 2)));
  
    int _numberOfSegments = ceil(trackLength / _segmentLength );
    double dEmean = (dd4hep::keV * _energyLoss * trackLength) / ((double)_numberOfSegments);

    // TODO _segmentLength may be different from segmentLength, is it ok?
    double segmentLength = trackLength / ((double)_numberOfSegments);
    double _segmentDepth = _sensor.GetThickness() / ((double)_numberOfSegments);

    double z = -_sensor.GetHalfThickness() - 0.5 * _segmentDepth;

    double eSum = 0.0;
    vector<TimedSignalPoint> signal_buffer;

    for (int i = 0; i < _numberOfSegments; ++i)
    {
        // ionization point
        z += _segmentDepth;
        double x = pos[0] + tanx * (z - pos[2]);
        double y = pos[1] + tany * (z - pos[2]);
        // momentum in MeV/c, mass in MeV, tmax (delta cut) in MeV, 
        // length in mm, meanLoss eloss in MeV.
        double eloss = _fluctuate->SampleFluctuations(particleMomentum * dd4hep::keV / dd4hep::MeV,
                                                      particleMass * dd4hep::keV / dd4hep::MeV,
                                                      _cutOnDeltaRays,
                                                      segmentLength,
                                                      dEmean / dd4hep::MeV) * dd4hep::MeV;

        double DistanceToPlane = _sensor.GetHalfThickness() - z;
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
        //double SigmaDiff = sqrt(DriftLength / _sensor.GetThickness()) * _diffusionCoefficient;
        double SigmaDiff = DistanceToPlane * _diffusionCoefficient;

        double SigmaX = SigmaDiff * sqrt(1.0 + pow(_tanLorentzAngleX, 2));
        double SigmaY = SigmaDiff * sqrt(1.0 + pow(_tanLorentzAngleY, 2));

        // energy is in keV       
        double charge = (eloss / dd4hep::keV) * _electronsPerKeV;

        signal_buffer.push_back({
            xOnPlane,
            yOnPlane,
            SigmaX,
            SigmaY,
            charge,
            hit
        });

        eSum += eloss;
    }

    double hEdep = hit->getEDep() / dd4hep::GeV;
    // deltaEne is a charge??
    const double thr = _deltaEne / _electronsPerKeV * dd4hep::keV;
    while (hEdep > eSum + thr)
    {
        // Add additional charge sampled from an 1 / n^2 distribution.
        const double       q = randomTail( thr, hEdep - eSum );
        const unsigned int h = floor(RandFlat::shoot(0.0, (double)_numberOfSegments));
        signal_buffer[h].charge += q * _electronsPerKeV / dd4hep::keV;
        eSum += q;
    }

    if (streamlog::out.write<streamlog::DEBUG5>() || streamlog::out.write<streamlog::DEBUG6>())
#pragma omp critical
    {
        streamlog::out() << "Ionization Points:" << std::endl;
        if (streamlog::out.write<streamlog::DEBUG5>())
        {
            streamlog::out() << "Number of ionization points: " << _numberOfSegments
                             << ", G4 EDep = "  << hEdep << std::endl
                             << "Padding each segment charge (1/n^2 pdf) until total below "
                             << _deltaEne << "e- threshold. New total energy: "
                             << eSum << std::endl;
        }
        else
        {
            streamlog::out() <<  "Track path length: " << trackLength
                             << ", calculated dEmean * N_segment = " << dEmean
                             << " * " << _numberOfSegments << " = "
                             << dEmean*_numberOfSegments << std::endl;
        }
    }

    for(auto spoint : signal_buffer)
    {
        signals.push_back(spoint);
    }
}

//=============================================================================
// Sample charge from 1 / n^2 distribution.
//=============================================================================
double DetElemSlidingWindow::randomTail( const double qmin, const double qmax )
{
    const double offset = 1. / qmax;
    const double range  = ( 1. / qmin ) - offset;
    const double u      = offset + RandFlat::shoot() * range;
    return 1. / u;
}


