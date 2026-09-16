#include "MuonCVXDDigitiser.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <sstream>

#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/Detector.h"
#include "DDRec/Surface.h"
#include "gsl/gsl_math.h"
#include "gsl/gsl_sf_erf.h"

#include "edm4hep/MCParticle.h"
#include "edm4hep/MutableSimTrackerHit.h"
#include "edm4hep/MutableTrackerHitPlane.h"

#include "EventRandom.h"

using dd4hep::rec::SurfaceMap;
using dd4hep::rec::ISurface;
using dd4hep::rec::Vector2D;
using dd4hep::rec::Vector3D;

DECLARE_COMPONENT(MuonCVXDDigitiser)

MuonCVXDDigitiser::MuonCVXDDigitiser(const std::string& name, ISvcLocator* svcLoc) :
    MultiTransformer(name, svcLoc,
                     {KeyValues("CollectionName", {"VertexBarrelCollection"}),
                      KeyValues("EventHeader", {"EventHeader"})},
                     {KeyValues("SimHitLocCollectionName", {"VertexBarrelPixels"}),
                      KeyValues("OutputCollectionName", {"VTXTrackerHits"}),
                      KeyValues("RelationColName", {"VTXTrackerHitRelations"}),
                      KeyValues("RawHitsLinkColName", {"VTXRawHitRelations"})})
{}

StatusCode MuonCVXDDigitiser::initialize()
{
    StatusCode sc = MultiTransformer::initialize();
    if (sc.isFailure()) return sc;

    m_uidSvc = service<IUniqueIDGenSvc>("UniqueIDGenSvc", true);
    if (!m_uidSvc) {
        error() << "Unable to get UniqueIDGenSvc" << endmsg;
        return StatusCode::FAILURE;
    }
    m_geoSvc = serviceLocator()->service(m_geoSvcName);
    if (!m_geoSvc) {
        error() << "Unable to retrieve the GeoSvc" << endmsg;
        return StatusCode::FAILURE;
    }

    try {
        // The pixel matrix is the one of the sensor the hit is on: z segmentation of the ladders does not matter
        m_geo = loadLayerGeometry(*m_geoSvc->getDetector(), m_subDetName, -1, m_layerIDs);
        m_cellID = std::make_unique<TrackerCellID>(m_geoSvc->constantAsString(m_encodingStringVariable));
    } catch (const std::exception& ex) {
        error() << ex.what() << endmsg;
        return StatusCode::FAILURE;
    }

    // Determine if vertex, inner tracker, or outer tracker
    if (!m_geo.type.isVertex && !m_geo.type.isInnerTracker && !m_geo.type.isOuterTracker) {
        error() << "Could not determine sub-detector type for: " << m_subDetName.value() << endmsg;
        return StatusCode::FAILURE;
    }
    if (m_layerIDs.empty()) {
        error() << "LayerIDs is empty: no SimTrackerHit could be assigned to a layer of " << m_subDetName.value() << endmsg;
        return StatusCode::FAILURE;
    }
    if ((int)m_layerIDs.size() != m_geo.numberOfLayers) {
        warning() << "LayerIDs has " << m_layerIDs.size() << " entries but " << m_subDetName.value() << " has "
                  << m_geo.numberOfLayers << " layers" << endmsg;
    }
    if (m_TimeDigitizeBinning != 0) {
        warning() << "Invalid setting for pixel time digitization binning (" << m_TimeDigitizeBinning.value()
                  << "): original times will be retained" << endmsg;
    }

    PrintGeometryInfo();
    FillChargeBins();

    return StatusCode::SUCCESS;
}

void MuonCVXDDigitiser::FillChargeBins()
{
    // Bins for charge discretization
    // FIXME: Will move to assign more dynamically
    std::vector<double> genericBins;
    if (m_ChargeDigitizeNumBits == 3) genericBins = {500, 786, 1100, 1451, 1854, 2390, 3326, 31973};
    if (m_ChargeDigitizeNumBits == 4) genericBins = {500, 657, 862, 1132, 1487, 1952, 2563, 3366, 4420, 5804, 7621, 10008, 13142, 17257, 22660, 29756};
    if (m_ChargeDigitizeNumBits == 5) genericBins = {500, 573, 633, 698, 757, 821, 890, 963, 1032, 1104, 1179, 1260, 1337, 1421, 1505, 1600, 1685, 1777, 1875, 1982, 2097, 2220, 2352, 2511, 2679, 2866, 3107, 3429, 3880, 4618, 6287, 16039};
    if (m_ChargeDigitizeNumBits == 6) genericBins = {500, 542, 572, 601, 629, 661, 692, 721, 750, 779, 812, 842, 877, 913, 946, 981, 1016, 1051, 1087, 1121, 1161, 1196, 1237, 1275, 1313, 1350, 1391, 1431, 1468, 1514, 1560, 1606, 1646, 1687, 1733, 1777, 1821, 1872, 1920, 1976, 2036, 2091, 2145, 2213, 2272, 2337, 2411, 2488, 2573, 2651, 2739, 2834, 2938, 3053, 3194, 3356, 3532, 3764, 4034, 4379, 4907, 5698, 6957, 9636};
    if (m_ChargeDigitizeNumBits == 8) genericBins = {500, 511, 523, 533, 542, 550, 556, 564, 570, 577, 585, 592, 598, 603, 610, 617, 624, 630, 638, 646, 654, 661, 668, 676, 684, 691, 699, 705, 712, 719, 724, 731, 738, 745, 752, 760, 767, 772, 780, 787, 795, 802, 810, 818, 826, 832, 839, 847, 856, 865, 874, 881, 889, 898, 906, 916, 924, 930, 938, 945, 955, 965, 971, 978, 986, 995, 1004, 1012, 1019, 1027, 1036, 1044, 1053, 1062, 1071, 1079, 1088, 1096, 1104, 1112, 1121, 1131, 1139, 1149, 1158, 1168, 1175, 1184, 1193, 1203, 1211, 1221, 1233, 1241, 1249, 1259, 1268, 1277, 1286, 1294, 1303, 1313, 1321, 1330, 1338, 1348, 1357, 1368, 1378, 1387, 1395, 1406, 1417, 1426, 1434, 1445, 1452, 1460, 1470, 1480, 1492, 1503, 1514, 1525, 1536, 1550, 1560, 1570, 1580, 1592, 1604, 1614, 1623, 1634, 1644, 1653, 1662, 1673, 1684, 1695, 1707, 1717, 1727, 1737, 1747, 1759, 1769, 1780, 1790, 1800, 1812, 1823, 1835, 1846, 1860, 1873, 1885, 1897, 1907, 1918, 1931, 1943, 1958, 1971, 1987, 2000, 2014, 2026, 2041, 2056, 2068, 2080, 2095, 2108, 2119, 2131, 2147, 2162, 2180, 2195, 2213, 2224, 2238, 2256, 2269, 2284, 2300, 2314, 2332, 2351, 2366, 2383, 2401, 2421, 2440, 2458, 2475, 2496, 2519, 2538, 2559, 2581, 2601, 2618, 2636, 2658, 2681, 2703, 2722, 2742, 2767, 2791, 2811, 2836, 2857, 2884, 2913, 2938, 2967, 2995, 3023, 3052, 3086, 3119, 3153, 3188, 3221, 3270, 3304, 3342, 3390, 3428, 3473, 3515, 3556, 3611, 3691, 3742, 3801, 3857, 3928, 3999, 4069, 4141, 4220, 4325, 4417, 4518, 4655, 4789, 4965, 5141, 5359, 5548, 5770, 6017, 6311, 6584, 7024, 7492, 8060, 8740, 9738, 11450, 14878, 23973};

    // Dedicated 4-bit tables for the vertex detector (barrel and endcap), keyed on the
    // sensor thickness in microns. Keyed on a rounded integer because layerThickness is
    // float and would never compare equal to a double literal.
    const std::map<int, std::vector<double>> vertexBins4bit = {
        {50, {500, 657, 862, 1132, 1487, 1952, 2563, 3366, 4420, 5804, 7621, 10008, 13142, 17257, 22660, 29756}},
        {75, {500, 675, 910, 1228, 1656, 2235, 3015, 4067, 5487, 7403, 9987, 13473, 18177, 24523, 33084, 44634}},
        {100, {500, 688, 946, 1300, 1788, 2460, 3382, 4652, 6397, 8797, 12098, 16638, 22881, 31467, 43274, 59512}},
        {200, {500, 720, 1037, 1494, 2152, 3099, 4463, 6428, 9258, 13334, 19205, 27660, 39838, 57378, 82640, 119024}},
        {400, {500, 754, 1138, 1716, 2588, 3904, 5889, 8883, 13399, 20211, 30486, 45985, 69363, 104626, 157816, 238048}},
    };

    // One table per layer, chosen from the sensor thickness of the layer
    m_DigitizedBins.assign(m_geo.numberOfLayers, genericBins);
    for (int i = 0; i < m_geo.numberOfLayers; ++i)
    {
        int thickness_um = (int)std::lround(m_geo.layerThickness[i] * 1000.);
        if (m_ChargeDigitizeNumBits == 4 && m_geo.type.isVertex)
        {
            auto it = vertexBins4bit.find(thickness_um);
            if (it != vertexBins4bit.end()) m_DigitizedBins[i] = it->second;
        }
        if (m_DigitizedBins[i].empty())
        {
            warning() << "No charge digitization bins defined for " << m_ChargeDigitizeNumBits.value()
                      << " bits; variable binning is unusable for layer " << i << endmsg;
            continue;
        }
        debug() << "Layer " << i << ": sensor thickness " << thickness_um
                << " um, " << m_DigitizedBins[i].size() << " charge bins, first/last = "
                << m_DigitizedBins[i].front() << "/" << m_DigitizedBins[i].back() << endmsg;
    }

    // shift digitized bins for inner and outer tracker by factor of 2
    // this adjusts for the fact that the resolution is 2x worse for inner and outer tracker
    if (!m_geo.type.isVertex) {
        debug() << "Subdetector is: " << m_subDetName.value() << endmsg;
        float shift = 500.; // first bin
        float scalefactor = 2.;
        for (auto& bins : m_DigitizedBins) {
            for (int i = 0; i < (int)bins.size(); i++) {
                bins[i] = (bins[i] - bins[0]) * scalefactor + shift;
            }
        }
    }
}

std::tuple<edm4hep::SimTrackerHitCollection,
           edm4hep::TrackerHitPlaneCollection,
           edm4hep::TrackerHitSimTrackerHitLinkCollection,
           edm4hep::TrackerHitSimTrackerHitLinkCollection>
MuonCVXDDigitiser::operator()(const edm4hep::SimTrackerHitCollection& STHcol,
                              const edm4hep::EventHeaderCollection& headers) const
{
    //SP. few TODO items:
    // - include noisy pixels (calculate rate from gaussian with unit sigma integral x > _electronicNoise / _threshold )
    // - change logic in creating pixels from all SimTrkHits, then cluster them (incl. timing info)
    // - include threshold dispersion effects
    // - add digi parametrization for time measurement
    // - change position determination of cluster to analog cluster (w-avg of corner hits)
    edm4hep::SimTrackerHitCollection STHLocCol;
    edm4hep::TrackerHitPlaneCollection THcol;
    edm4hep::TrackerHitSimTrackerHitLinkCollection relCol;
    edm4hep::TrackerHitSimTrackerHitLinkCollection rawHitsCol;

    const auto seed = m_uidSvc->getUniqueID(headers, name());
    debug() << "Using seed " << seed << " for event " << headers[0].getEventNumber() << " and run "
            << headers[0].getRunNumber() << endmsg;
    EventRandom random(seed);

    const bool debugOn = msgLevel(MSG::DEBUG);
    const bool verboseOn = msgLevel(MSG::VERBOSE);

    HitState state;
    PixelHitMap pixels;

    int nSimHits = STHcol.size();
    debug() << "Processing collection with " << nSimHits << " hits ... " << endmsg;

    for (int i = 0; i < nSimHits; ++i)
    {
        const edm4hep::SimTrackerHit simTrkHit = STHcol[i];
        const auto cellid = simTrkHit.getCellID();
        const auto& simPos = simTrkHit.getPosition();
        const auto& simMom = simTrkHit.getMomentum();
        // use CellID to set layer and ladder numbers
        state.currentLayer = layerMapping(m_cellID->layer(cellid));
        state.currentLadder = m_cellID->module(cellid);
        if (debugOn)
        {
            float sim_r = std::sqrt(simPos.x * simPos.x + simPos.y * simPos.y);
            float sim_phi = std::atan(simPos.y / simPos.x);
            float sim_theta = simPos.z == 0 ? 3.1416 / 2 : std::atan(sim_r / simPos.z);
            debug() << "Processing simHit #" << i << ", from layer=" << state.currentLayer << ", module=" << state.currentLadder << "\n"
                    << "- EDep = " << simTrkHit.getEDep() * dd4hep::GeV / dd4hep::keV << " keV, path length = " << simTrkHit.getPathLength() * 1000. << " um\n"
                    << "- Position (mm) x,y,z,t = " << simPos.x << ", " << simPos.y << ", " << simPos.z << ", " << simTrkHit.getTime() << "\n"
                    << "- Position r(mm),phi,theta = " << sim_r << ", " << sim_phi << ", " << sim_theta << "\n"
                    << "- MC particle pdg = ";
            const edm4hep::MCParticle mcp = simTrkHit.getParticle();
            if (mcp.isAvailable()) {
                debug() << mcp.getPDG() << "\n";
                double deltaZ = std::fabs(simPos.z - mcp.getVertex().z);
                double deltaX = simPos.x - mcp.getVertex().x;
                double deltaY = simPos.y - mcp.getVertex().y;
                double deltaR = std::sqrt(deltaX * deltaX + deltaY * deltaY);
                double incidentTheta = std::atan2(deltaR, deltaZ);
                double p = std::sqrt(std::pow(simMom.x, 2) + std::pow(simMom.y, 2) + std::pow(simMom.z, 2));
                double mass = mcp.getMass();
                debug() << "delta r: " << deltaR << "\n"
                        << "delta z: " << deltaZ << "\n"
                        << "incident theta: " << incidentTheta << " radians or " << incidentTheta * (180 / M_PI) << " degrees\n"
                        << "beta: " << p / std::sqrt(p * p + mass * mass) << "\n";
            } else {
                debug() << " N.A.\n";
            }
            debug() << "- MC particle p (GeV) = " << std::sqrt(simMom.x * simMom.x + simMom.y * simMom.y + simMom.z * simMom.z) << "\n"
                    << "- isSecondary = " << simTrkHit.isProducedBySecondary() << ", isOverlay = " << simTrkHit.isOverlay() << "\n"
                    << "- Quality = " << simTrkHit.getQuality() << endmsg;
        }

        if (state.currentLayer == -1)
            continue;
        ProduceIonisationPoints(simTrkHit, state, random);
        if (state.currentLayer == -1) {
            ++m_nOffSurface;
            continue;
        }
        ProduceSignalPoints(state);
        ProduceHits(pixels, simTrkHit, state);
        if (m_PoissonSmearing) PoissonSmearer(pixels, random);
        if (m_electronicEffects) GainSmearer(pixels, random);
        ApplyThreshold(pixels, random);
        if (m_DigitizeCharge) ChargeDigitizer(pixels, state);
        if (m_timeSmearingModel != 0) TimeSmearer(pixels, state, random);
        if (m_DigitizeTime) TimeDigitizer(pixels);

        //**************************************************************************
        // Create reconstructed cluster object
        //**************************************************************************
        ClusterHit cluster;
        if (!ReconstructTrackerHit(pixels, state, cluster))
        {
            debug() << "Skip hit" << endmsg;
            ++m_nNoCharge;
            continue;
        }

        auto recoHit = THcol.create();
        recoHit.setEDep(cluster.eDep);
        recoHit.setDu(m_pixelSizeX / std::sqrt(12));
        recoHit.setDv(m_pixelSizeY / std::sqrt(12));
        recoHit.setTime(cluster.time);
        // hit's layer/ladder/petal position does not change
        recoHit.setCellID(cellid);

        double xLab[3];
        TransformToLab(cellid, cluster.position, xLab);
        recoHit.setPosition({xLab[0], xLab[1], xLab[2]});

        const ISurface* surf = m_geo.surfaceMap->find(cellid)->second;
        Vector3D u = surf->u();
        Vector3D v = surf->v();
        //TODO HACK: Store incidence angle of particle instead!
        recoHit.setU({float(u.theta()), float(u.phi())});
        recoHit.setV({float(v.theta()), float(v.phi())});

        if (debugOn)
        {
            // Local position of the SimTrackerHit, only needed to check the reconstruction
            double localPos[3];
            double localDir[3];
            HitState checkState = state;
            FindLocalPosition(simTrkHit, localPos, localDir, checkState);
            float incidentPhi = std::atan(localDir[0] / localDir[2]);
            float incidentTheta = std::atan(localDir[1] / localDir[2]);
            const auto& recoPos = recoHit.getPosition();
            debug() << "- TRUE GLOBAL position (mm) x,y,z,t = " << simPos.x << ", " << simPos.y << ", " << simPos.z << ", " << simTrkHit.getTime() << "\n"
                    << "- TRUE LOCAL position (localPos) (mm) x,y,z,t = " << localPos[0] << ", " << localPos[1] << ", " << localPos[2] << "\n"
                    << "- RECO LOCAL position (mm) x,y,z,t = " << cluster.position[0] << ", " << cluster.position[1] << ", " << cluster.position[2] << "\n"
                    << "- RECO GLOBAL position (mm) x,y,z,t = " << recoPos.x << ", " << recoPos.y << ", " << recoPos.z << "\n"
                    << "Reconstructed pixel cluster:\n"
                    << "- local position (x,y) = " << localPos[0] << "(Idx: " << localPos[0] / m_pixelSizeX << "), "
                    << localPos[1] << "(Idy: " << localPos[1] / m_pixelSizeY << ")\n"
                    << "(reco local) - (true local) (x,y,z): " << localPos[0] - state.currentLocalPosition[0] << ", "
                    << localPos[1] - state.currentLocalPosition[1] << ", " << localPos[2] - state.currentLocalPosition[2] << "\n"
                    << "- global position (x,y,z, t) = " << recoPos.x << ", " << recoPos.y << ", " << recoPos.z << ", " << recoHit.getTime() << "\n"
                    << "- (reco global (x,y,z,t)) - (true global) = " << recoPos.x - simPos.x << ", " << recoPos.y - simPos.y << ", "
                    << recoPos.z - simPos.z << ", " << recoHit.getTime() - simTrkHit.getTime() << "\n"
                    << "- charge = " << recoHit.getEDep() << "(True: " << simTrkHit.getEDep() << ")\n"
                    << "- incidence angles: theta = " << incidentTheta << ", phi = " << incidentPhi << endmsg;
        }

        //**************************************************************************
        // Set Relation to SimTrackerHit
        //**************************************************************************
        auto rel = relCol.create();
        rel.setFrom(recoHit);
        rel.setTo(simTrkHit);
        rel.setWeight(1.0);

        if (m_produceFullPattern)
        {
            // Store all the fired points
            std::vector<edm4hep::MutableSimTrackerHit> rawHits;
            for (const auto& [pixelID, sth] : pixels)
            {
                float charge = sth.charge;
                //store hits that are above threshold. In case of _ChargeDiscretization, just check for a small non-zero value
                if ((m_DigitizeCharge and (charge > 1.0)) or (charge > m_threshold))
                {
                    auto newsth = STHLocCol.create();
                    // hit's layer/ladder position is the same for all fired points
                    newsth.setCellID(cellid);
                    //Store local position in units of pixels instead
                    newsth.setPosition({sth.x / m_pixelSizeX, sth.y / m_pixelSizeY, 0.});
                    newsth.setEDep(charge); // in unit of electrons
                    newsth.setTime(sth.time);
                    newsth.setPathLength(simTrkHit.getPathLength());
                    newsth.setParticle(simTrkHit.getParticle());
                    newsth.setMomentum(simMom);
                    newsth.setProducedBySecondary(simTrkHit.isProducedBySecondary());
                    newsth.setOverlay(simTrkHit.isOverlay());
                    rawHits.push_back(newsth);
                }
            }
            for (const auto& rawHit : rawHits)
            {
                auto rawLink = rawHitsCol.create();
                rawLink.setFrom(recoHit);
                rawLink.setTo(rawHit);
                rawLink.setWeight(1. / rawHits.size());
            }
            if (verboseOn)
            {
                verbose() << "- number of pixels: " << rawHits.size() << "\n"
                          << "- List of constituents (pixels/strips):\n";
                for (size_t iH = 0; iH < rawHits.size(); ++iH)
                    verbose() << "  - " << iH << ": Edep (e-) = " << rawHits[iH].getEDep() << ", t (ns) =" << rawHits[iH].getTime() << "\n";
                verbose() << "--------------------------------" << endmsg;
            }
        }
    }
    debug() << "Number of produced hits: " << THcol.size() << endmsg;

    m_nSimHits += STHcol.size();
    m_nRecoHits += THcol.size();

    return std::make_tuple(std::move(STHLocCol), std::move(THcol), std::move(relCol), std::move(rawHitsCol));
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
void MuonCVXDDigitiser::FindLocalPosition(const edm4hep::SimTrackerHit& hit,
                                          double* localPosition,
                                          double* localDirection,
                                          HitState& state) const
{
    // Use SurfaceManager to calculate local coordinates
    const auto cellID = hit.getCellID();
    verbose() << "Cell ID of Sim Hit: " << cellID << endmsg;
    SurfaceMap::const_iterator sI = m_geo.surfaceMap->find(cellID);
    if (sI == m_geo.surfaceMap->end()) {
        debug() << "  no surface found for cell ID " << cellID << endmsg;
        state.currentLayer = -1;
        return;
    }
    const ISurface* surf = sI->second;
    // The pixel matrix covers the sensor the hit is on
    const SensorExtent* extent = m_geo.sensorExtent(cellID);
    state.sensorHalfLengthU = extent->halfLengthU;
    state.sensorHalfLengthV = extent->halfLengthV;
    Vector3D oldPos(hit.getPosition().x, hit.getPosition().y, hit.getPosition().z);
    // We need it?
    if (!surf->insideBounds(dd4hep::mm * oldPos)) {
        if (msgLevel(MSG::VERBOSE)) {
            std::ostringstream surfDesc;
            surfDesc << *surf;
            verbose() << "  hit at " << oldPos
                      << " is not on surface "
                      << surfDesc.str()
                      << " distance: " << surf->distance(dd4hep::mm * oldPos)
                      << endmsg;
        }
        state.currentLayer = -1;
        return;
    }

    Vector2D lv = surf->globalToLocal(dd4hep::mm * oldPos);
    // Store local position in mm
    localPosition[0] = lv[0] / dd4hep::mm;
    localPosition[1] = lv[1] / dd4hep::mm;
    // Add also z ccordinate
    Vector3D origin(surf->origin()[0], surf->origin()[1], surf->origin()[2]);
    localPosition[2] = (dd4hep::mm * oldPos - dd4hep::cm * origin).dot(surf->normal()) / dd4hep::mm;
    // Prefer the momentum recorded at the hit: that is the direction inside the sensor, with
    // everything upstream (bending, scattering, energy loss) already folded in. The MC
    // particle's momentum is taken at its production vertex, which in a solenoid can point
    // somewhere else entirely -- in a muon-gun sample the two differ by 3.5 deg at the median
    // but by more than 65 deg for a quarter of the hits. Fall back to the MC particle only
    // when the hit carries no momentum, as not every producer fills it.
    double Momentum[3];
    const edm4hep::MCParticle mcp = hit.getParticle();
    const auto& hitMomentum = hit.getMomentum();
    bool useHitMomentum = (hitMomentum.x != 0.f || hitMomentum.y != 0.f || hitMomentum.z != 0.f);
    const double hitMom[3] = {hitMomentum.x, hitMomentum.y, hitMomentum.z};
    for (int j = 0; j < 3; ++j) {
        if (useHitMomentum) {
            Momentum[j] = hitMom[j] * dd4hep::GeV;
        } else if (mcp.isAvailable()) {
            const auto& mcpMomentum = mcp.getMomentum();
            const double mcpMom[3] = {mcpMomentum.x, mcpMomentum.y, mcpMomentum.z};
            Momentum[j] = mcpMom[j] * dd4hep::GeV;
        } else {
            Momentum[j] = 0.;
        }
    }
    if (!useHitMomentum)
        debug() << "Hit carries no momentum, falling back to the MC particle" << endmsg;
    // as default put electron's mass
    state.currentParticleMass = 0.510e-3 * dd4hep::GeV;
    if (mcp.isAvailable())
        state.currentParticleMass = std::max(mcp.getMass() * dd4hep::GeV, state.currentParticleMass);
    state.currentParticleMomentum = sqrt(pow(Momentum[0], 2) + pow(Momentum[1], 2)
                                         + pow(Momentum[2], 2));

    localDirection[0] = Momentum * surf->u();
    localDirection[1] = Momentum * surf->v();
    localDirection[2] = Momentum * surf->normal();
    if (m_geo.type.isBarrel) {
        state.currentPhi = state.currentLadder * 2.0 * m_geo.layerHalfPhi[state.currentLayer] + m_geo.layerPhiOffset[state.currentLayer];
    }
}

void MuonCVXDDigitiser::ProduceIonisationPoints(const edm4hep::SimTrackerHit& hit, HitState& state, EventRandom& random) const
{
    verbose() << "Creating Ionization Points" << endmsg;
    double pos[3] = {0, 0, 0};
    double dir[3] = {0, 0, 0};
    double entry[3];
    double exit[3];
    // hit and pos are in mm
    FindLocalPosition(hit, pos, dir, state);
    if (state.currentLayer == -1)
        return;

    const int layer = state.currentLayer;
    double origPos[3] = {pos[0], pos[1], pos[2]};
    entry[2] = -m_geo.layerHalfThickness[layer];
    exit[2] = m_geo.layerHalfThickness[layer];
    // entry points: hit position is in middle of layer. ex: entry_x = x - (z distance to bottom of layer) * px/pz
    for (int i = 0; i < 2; ++i) {
        entry[i] = origPos[i] + dir[i] * (entry[2] - origPos[2]) / dir[2];
        exit[i] = origPos[i] + dir[i] * (exit[2] - origPos[2]) / dir[2];
    }

    // PDG multiple-scattering charge number z: |q| in units of e. Fall back to unit charge
    // when the hit carries no MC particle, matching the electron-mass fallback in
    // FindLocalPosition(). Neutrals are skipped: with z = 0 the log term is log(0) and
    // theta_0 would evaluate to 0 * -inf = NaN.
    const edm4hep::MCParticle ms_mcp = hit.getParticle();
    double q_charge = ms_mcp.isAvailable() ? std::fabs(ms_mcp.getCharge()) : 1.0;

    if (m_doMultipleScattering && q_charge == 0.)
        verbose() << "Neutral particle, skipping multiple scattering" << endmsg;

    if (m_doMultipleScattering && q_charge > 0.) {
        verbose() << "Applying multiple scattering formula" << endmsg;
        //**************************************************************************
        //Multiple scattering implementation based on PDG formula (https://pdg.lbl.gov/2020/reviews/rpp2020-rev-passage-particles-matter.pdf)
        //**************************************************************************
        const auto& hitMomentum = hit.getMomentum();
        double p = std::sqrt(pow(hitMomentum.x, 2) + pow(hitMomentum.y, 2) + pow(hitMomentum.z, 2)); //[GeV/c]
        double beta = p / std::sqrt(p * p + state.currentParticleMass * state.currentParticleMass);

        double x_0 = 93.7; // [mm] -> radiation length in silicon
        double sensorT = m_geo.layerThickness[layer]; // [mm] -> sensor thickness

        // normalize direction
        double mag = std::sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);
        dir[0] /= mag;
        dir[1] /= mag;
        dir[2] /= mag;
        for (int i = 0; i < 2; ++i) {
            entry[i] = pos[i] + dir[i] * (entry[2] - pos[2]) / dir[2];
        }

        double pathL_segment, theta_0, theta_plane_x, theta_plane_y, theta_out_x, theta_out_y;
        // Step through the sensor in a whole number of slices no thicker than
        // MSSliceThickness, so the steps cover exactly the sensor thickness. Accumulating a
        // floating-point z against sensorT instead runs a whole extra slice: layerThickness
        // is a float, so 50 um / 5 um evaluates to 10.000000149 rather than 10. The small
        // relative tolerance absorbs that while still rounding a genuine 10.4 up to 11.
        double sliceT = (m_msSliceThickness > 0.) ? m_msSliceThickness.value() : sensorT;
        int n_slices = std::max(1, (int)std::ceil(sensorT / sliceT * (1. - 1e-6)));
        double z_segment = sensorT / n_slices;
        for (int islice = 0; islice < n_slices; ++islice) {

            pathL_segment = z_segment / fabs(dir[2]); // path length for segment

            // --- Multiple scattering step ---
            // 0.0136 GeV is the PDG 13.6 MeV constant; natural units (c = 1) with p in GeV/c
            theta_0 = (0.0136 / (beta * p)) * q_charge * std::sqrt(pathL_segment / x_0)
                * (1 + 0.038 * std::log(pathL_segment * std::pow(q_charge, 2) / (x_0 * std::pow(beta, 2)))); //as defined in PDG

            theta_plane_x = random.gauss(0., theta_0);
            theta_plane_y = random.gauss(0., theta_0);
            theta_out_x = theta_plane_x + std::atan2(dir[0], dir[2]);
            theta_out_y = theta_plane_y + std::atan2(dir[1], dir[2]);

            //update dir vector:
            dir[0] = std::tan(theta_out_x);
            dir[1] = std::tan(theta_out_y);
            dir[2] = 1.0;

            // renormalize again
            double mag2 = std::sqrt(dir[0] * dir[0] + dir[1] * dir[1] + 1.0);
            dir[0] /= mag2;
            dir[1] /= mag2;
            dir[2] /= mag2;

            //update postion:
            pos[0] += dir[0] * pathL_segment;
            pos[1] += dir[1] * pathL_segment;
            pos[2] += dir[2] * pathL_segment;
        }
        //find final exit point
        for (int i = 0; i < 2; ++i) {
            exit[i] = pos[i] + dir[i] * (exit[2] - pos[2]) / dir[2];
        }
    }
    //end of multiple scattering implementation

    for (int i = 0; i < 3; ++i) {
        state.currentLocalPosition[i] = origPos[i];
        state.currentEntryPoint[i] = entry[i];
        state.currentExitPoint[i] = exit[i];
    }
    debug() << "local position: " << state.currentLocalPosition[0] << ", " << state.currentLocalPosition[1] << ", " << state.currentLocalPosition[2] << endmsg;
    // Local unit direction of travel through the sensor. After the multiple-scattering loop
    // above this is the deflected direction, so the re-simulation chain sees its own kinematics.
    double dirMag = std::sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);
    double u[3] = {dir[0] / dirMag, dir[1] / dirMag, dir[2] / dirMag};
    const double halfT = m_geo.layerHalfThickness[layer];
    const bool crossesPlane = std::fabs(u[2]) > 1e-12;

    // Straight-line crossing of the full sensor along that direction.
    double crossingLength = crossesPlane
                          ? m_geo.layerThickness[layer] / std::fabs(u[2])
                          : m_maxTrkLen.value();

    // Two mutually exclusive sources for the trail, selected by ResimulateIonisation:
    //  false (Geant4, default): the trail length is the path length Geant4 recorded in the
    //    sensitive volume and the total deposit is the Geant4 EDep.
    //  true (re-simulation): the trail is the full crossing of the sensor along the local
    //    direction and the energy comes from the EnergyLoss dE/dx parametrisation. Nothing
    //    but the direction is taken from Geant4, so the chain can be tuned on its own.
    const bool resimulate = m_resimulateIonisation;

    double pathLength;
    if (resimulate)
    {
        pathLength = crossingLength;
    }
    else
    {
        // True path length in the sensitive volume as recorded by Geant4: it already contains
        // the scattering and curvature inside the sensor, and it is shorter than a full
        // crossing for a particle that stopped, started or clipped a corner. Fall back to the
        // straight-line crossing when the producer did not fill it.
        pathLength = hit.getPathLength();
        if (!(pathLength > 0.))
        {
            pathLength = crossingLength;
            debug() << "Hit carries no path length, falling back to the sensor crossing" << endmsg;
        }
    }

    // trackLength is the physical distance travelled: it sets the segmentation and the
    // per-segment path length handed to the fluctuation model. Limited at 1cm.
    double trackLength = std::min(m_maxTrkLen.value(), pathLength);

    state.numberOfSegments = ceil(trackLength / m_segmentLength);

    // Energy Geant4 actually deposited in the sensor. In G4 mode it is the truth for the total
    // and anchors both dEmean and the 1/n^2 padding below; SampleFluctuations() still supplies
    // the segment-to-segment Landau structure around it. In re-simulation mode the total is
    // left to the parametrisation and the padding is skipped, so the hit stays independent of
    // the Geant4 deposit. The parametrisation is also the fallback for hits with no deposit.
    double hcharge = hit.getEDep() * dd4hep::GeV;
    const bool anchorToG4 = !resimulate && (hcharge > 0.);
    double dEmean = anchorToG4
                  ? hcharge / ((double)state.numberOfSegments)
                  : (dd4hep::keV * m_energyLoss * trackLength) / ((double)state.numberOfSegments);
    state.ionisationPoints.resize(state.numberOfSegments);
    verbose() << "Track path length: " << trackLength << ", calculated dEmean * N_segment = " << dEmean << " * "
              << state.numberOfSegments << " = " << dEmean * state.numberOfSegments << endmsg;
    state.eSum = 0.0;
    // TODO SegmentLength may be different from segmentLength, is it ok?
    double segmentLength = trackLength / ((double)state.numberOfSegments);

    // Place the trail as a straight segment along u, parametrised by the path coordinate s
    // measured from the hit position.
    double sLo, sHi;
    if (resimulate)
    {
        if (crossesPlane)
        {
            // The trail spans the whole sensor thickness, entry face to exit face.
            sLo = (-halfT - origPos[2]) / u[2];
            sHi = ( halfT - origPos[2]) / u[2];
            if (sLo > sHi) std::swap(sLo, sHi);
        }
        else
        {
            // Track running in the sensor plane: no crossing to span, keep the trail centred.
            sLo = -0.5 * trackLength;
            sHi = -sLo;
        }
    }
    else
    {
        // Centred on the hit position, its extent capped at the full-thickness crossing and
        // then clipped to the slab: a curling track reports an arc far longer than any straight
        // segment through the sensor, and ProduceSignalPoints() computes the drift distance as
        // (halfThickness - z), which goes negative if a point escapes.
        sLo = -0.5 * std::min(trackLength, crossingLength);
        sHi = -sLo;
        if (crossesPlane) {
            double sA = (-halfT - origPos[2]) / u[2];
            double sB = ( halfT - origPos[2]) / u[2];
            if (sA > sB) std::swap(sA, sB);
            sLo = std::max(sLo, sA);
            sHi = std::min(sHi, sB);
        }
    }
    if (!(sHi > sLo)) { sLo = 0.; sHi = 0.; }  // degenerate, put everything at the hit
    double geomStep = (sHi - sLo) / ((double)state.numberOfSegments);
    state.segmentDepth = geomStep * std::fabs(u[2]);

    debug() << "Number of ionization points: " << state.numberOfSegments
            << ", trail from " << (resimulate ? "re-simulation" : "Geant4")
            << ", G4 EDep = " << hcharge << endmsg;
    for (int i = 0; i < state.numberOfSegments; ++i)
    {
        double sPath = sLo + (i + 0.5) * geomStep;
        double x = origPos[0] + sPath * u[0];
        double y = origPos[1] + sPath * u[1];
        double z = origPos[2] + sPath * u[2];
        // momentum in MeV/c, mass in MeV, tmax (delta cut) in MeV,
        // length in mm, meanLoss eloss in MeV.
        // The delta-ray cut is copied for every call: SampleFluctuations() modifies it.
        double tmax = m_cutOnDeltaRays;
        double de = m_fluctuate.SampleFluctuations(random.engine(),
                                                   double(state.currentParticleMomentum / dd4hep::MeV),
                                                   double(state.currentParticleMass / dd4hep::MeV),
                                                   tmax,
                                                   segmentLength,
                                                   double(dEmean / dd4hep::MeV)) * dd4hep::MeV;
        state.eSum = state.eSum + de;
        IonisationPoint ipoint;
        ipoint.eloss = de;
        ipoint.x = x;
        ipoint.y = y;
        ipoint.z = z;
        state.ionisationPoints[i] = ipoint;
        verbose() << " " << i << ": z=" << z << ", eloss = " << de << "(total so far: " << state.eSum << "), x=" << x << ", y=" << y << endmsg;
    }

    // Top the sampled total up to the Geant4 deposit. Only meaningful when the hit is anchored
    // to Geant4: in re-simulation mode the total is whatever the parametrisation sampled.
    if (anchorToG4)
    {
        const double thr = m_deltaEne / m_electronsPerKeV * dd4hep::keV;
        while (hcharge > state.eSum + thr) {
            // Add additional charge sampled from an 1 / n^2 distribution.
            // Adjust charge to match expectations
            const double q = randomTail(thr, hcharge - state.eSum, random);
            const unsigned int h = floor(random.flat(0.0, (double)state.numberOfSegments));
            state.ionisationPoints[h].eloss += q;
            state.eSum += q;
        }
        debug() << "Padding each segment charge (1/n^2 pdf) until total below " << m_deltaEne.value()
                << "e- threshold. New total energy: " << state.eSum << endmsg;
    }
    if (msgLevel(MSG::VERBOSE)) {
        verbose() << "List of ionization points:\n";
        for (int i = 0; i < state.numberOfSegments; ++i) {
            verbose() << "- " << i << ": E=" << state.ionisationPoints[i].eloss
                      << ", x=" << state.ionisationPoints[i].x << ", y=" << state.ionisationPoints[i].y
                      << ", z=" << state.ionisationPoints[i].z << "\n";
        }
        verbose() << endmsg;
    }
}

void MuonCVXDDigitiser::ProduceSignalPoints(HitState& state) const
{
    state.signalPoints.resize(state.numberOfSegments);
    // run over ionisation points
    verbose() << "Creating signal points" << endmsg;
    for (int i = 0; i < state.numberOfSegments; ++i)
    {
        IonisationPoint ipoint = state.ionisationPoints[i]; // still local coords
        double z = ipoint.z;
        double x = ipoint.x;
        double y = ipoint.y;
        double DistanceToPlane = m_geo.layerHalfThickness[state.currentLayer] - z;
        double xOnPlane = x + m_tanLorentzAngleX * DistanceToPlane;
        double yOnPlane = y + m_tanLorentzAngleY * DistanceToPlane;
        // For diffusion-coeffieint calculation, see e.g. https://www.slac.stanford.edu/econf/C060717/papers/L008.PDF
        // or directly Eq. 13 of https://cds.cern.ch/record/2161627/files/ieee-tns-07272141.pdf
        // diffusionCoefficient = sqrt(2*D / mu / V), where
        //  - D = 12 cm^2/s // diffusion constant
        //  - mu = 450 cm^2/s/V // mobility
        //  - V = 10-30 V // expected depletion voltage
        //  => _diffusionCoefficient = 0.04-0.07
        // and diffusion sigma = _diffusionCoefficient * DistanceToPlane
        // e.g. fot 50um diffusion sigma = 2.1 - 3.7 um
        double SigmaDiff = DistanceToPlane * m_diffusionCoefficient;
        double SigmaX = SigmaDiff * sqrt(1.0 + pow(m_tanLorentzAngleX.value(), 2));
        double SigmaY = SigmaDiff * sqrt(1.0 + pow(m_tanLorentzAngleY.value(), 2));
        // energy is in keV
        double charge = (ipoint.eloss / dd4hep::keV) * m_electronsPerKeV;
        SignalPoint spoint;
        spoint.x = xOnPlane;
        spoint.y = yOnPlane;
        spoint.sigmaX = SigmaX;
        spoint.sigmaY = SigmaY;
        spoint.charge = charge; // electrons x keV
        state.signalPoints[i] = spoint;
        verbose() << "- " << i << ": charge=" << charge
                  << ", x=" << xOnPlane << "(delta=" << xOnPlane - x << ")"
                  << ", y=" << yOnPlane << "(delta=" << yOnPlane - y << ")"
                  << ", sigmaDiff=" << SigmaDiff
                  << ", sigmaX=" << SigmaX << ", sigmay=" << SigmaY << endmsg;
    }
}

void MuonCVXDDigitiser::ProduceHits(PixelHitMap& pixels, const edm4hep::SimTrackerHit& simHit, const HitState& state) const
{
    pixels.clear();
    const bool verboseOn = msgLevel(MSG::VERBOSE);
    verbose() << "Creating hits" << endmsg;
    const int pixelsInaColumn = GetPixelsInaColumn(state);
    const int pixelsInaRow = GetPixelsInaRow(state);
    for (int i = 0; i < state.numberOfSegments; ++i)
    {
        SignalPoint spoint = state.signalPoints[i];
        double xCentre = spoint.x;
        double yCentre = spoint.y;
        double sigmaX = spoint.sigmaX;
        double sigmaY = spoint.sigmaY;
        double xLo = spoint.x - 3 * spoint.sigmaX;
        double xUp = spoint.x + 3 * spoint.sigmaX;
        double yLo = spoint.y - 3 * spoint.sigmaY;
        double yUp = spoint.y + 3 * spoint.sigmaY;

        int ixLo, ixUp, iyLo, iyUp;
        TransformXYToCellID(xLo, yLo, ixLo, iyLo, state);
        TransformXYToCellID(xUp, yUp, ixUp, iyUp, state);
        verbose() << i << ": Pixel idx boundaries: ixLo=" << ixLo << ", iyLo=" << iyLo
                  << ", ixUp=" << ixUp << ", iyUp=" << iyUp << endmsg;
        for (int ix = ixLo; ix < ixUp + 1; ++ix)
        {
            if ((ix < 0) or (ix >= pixelsInaColumn)) {
                if (verboseOn) verbose() << "Pixels in a column: " << pixelsInaColumn << ", skipping pixels with ix =" << ix << endmsg;
                continue;
            }
            for (int iy = iyLo; iy < iyUp + 1; ++iy)
            {
                if ((iy < 0) or (iy >= pixelsInaRow)) {
                    if (verboseOn) verbose() << "Pixels in a row: " << pixelsInaRow << ", skipping pixels with iy =" << iy << endmsg;
                    continue;
                }
                double xCurrent, yCurrent;
                TransformCellIDToXY(ix, iy, xCurrent, yCurrent, state);

                gsl_sf_result result;
                /*int status = */gsl_sf_erf_Q_e((xCurrent - 0.5 * m_pixelSizeX - xCentre) / sigmaX, &result);
                double LowerBound = 1 - result.val;
                /*status = */gsl_sf_erf_Q_e((xCurrent + 0.5 * m_pixelSizeX - xCentre) / sigmaX, &result);
                double UpperBound = 1 - result.val;
                double integralX = UpperBound - LowerBound;
                /*status = */gsl_sf_erf_Q_e((yCurrent - 0.5 * m_pixelSizeY - yCentre) / sigmaY, &result);
                LowerBound = 1 - result.val;
                /*status = */gsl_sf_erf_Q_e((yCurrent + 0.5 * m_pixelSizeY - yCentre) / sigmaY, &result);
                UpperBound = 1 - result.val;
                double integralY = UpperBound - LowerBound;
                float totCharge = float(spoint.charge * integralX * integralY);
                int pixelID = pixelsInaRow * ix + iy;

                auto item = pixels.find(pixelID);
                if (item == pixels.end())
                {
                    // still in local coordinates, with the usual true timing as starting point
                    pixels.emplace(pixelID, PixelHit{xCurrent, yCurrent, m_geo.layerHalfThickness[state.currentLayer],
                                                     totCharge, simHit.getTime()});
                    if (verboseOn) verbose() << "Created new pixel hit at idx=" << ix << ", idy=" << iy << ", charge=" << totCharge << endmsg;
                }
                else
                {
                    item->second.charge += totCharge;
                    //TODO: handle multiple times. For now not needed since all deposits arrive at the same true time.
                    if (verboseOn) verbose() << "Updating pixel hit at idx=" << ix << ", idy=" << iy << ", total charge=" << item->second.charge
                              << "(delta = " << totCharge << ")" << endmsg;
                }
            }
        }
    }
    if (msgLevel(MSG::VERBOSE)) {
        verbose() << "List of pixel hits created:\n"; // still in local coords
        int idx = 0;
        for (const auto& [pixelID, hit] : pixels)
            verbose() << idx++ << ": x=" << hit.x << ", y=" << hit.y << ", z=" << hit.z << ", EDep = " << hit.charge << "\n";
        verbose() << endmsg;
    }
}

/**
 * Function that fluctuates charge (in units of electrons)
 * deposited on the fired pixels according to the Poisson
 * distribution...
 */
void MuonCVXDDigitiser::PoissonSmearer(PixelHitMap& pixels, EventRandom& random) const
{
    verbose() << "Adding Poisson smear to charge" << endmsg;
    for (auto& [pixelID, hit] : pixels)
    {
        float charge = hit.charge;
        float rng;
        if (charge > 1e+03) // assume Gaussian
        {
            rng = float(random.gauss(charge, sqrt(charge)));
        }
        else // assume Poisson
        {
            rng = float(random.poisson(charge));
        }
        hit.charge = rng;
        verbose() << pixelID << ": x=" << hit.x << ", y=" << hit.y << ", z=" << hit.z
                  << ", charge = " << rng << "(delta = " << charge - rng << ")" << endmsg;
    }
}

/**
 * Simulation of electronic noise.
 */
void MuonCVXDDigitiser::GainSmearer(PixelHitMap& pixels, EventRandom& random) const
{
    verbose() << "Adding FE noise smear to charge" << endmsg;
    for (auto& [pixelID, hit] : pixels)
    {
        double Noise = random.gauss(0., m_electronicNoise);
        hit.charge = hit.charge + float(Noise);
        verbose() << pixelID << ": x=" << hit.x << ", y=" << hit.y << ", z=" << hit.z
                  << ", charge = " << hit.charge << "(delta = " << Noise << ")" << endmsg;
    }
}

/**
 * Apply threshold.
 * Sets the charge to 0 if less than the threshold
 * Smears the threshold by a Gaussian if sigma > 0
 */
void MuonCVXDDigitiser::ApplyThreshold(PixelHitMap& pixels, EventRandom& random) const
{
    verbose() << "Applying threshold" << endmsg;
    for (auto& [pixelID, hit] : pixels)
    {
        double smear = 0;
        float origCharge = hit.charge;
        // Each pixel has its own threshold, smeared around the nominal one
        if (m_thresholdSmearSigma > 0) smear = random.gauss(0., m_thresholdSmearSigma);
        float actualThreshold = m_threshold + smear;
        if (hit.charge <= actualThreshold) hit.charge = 0.0;

        verbose() << pixelID << ": x=" << hit.x << ", y=" << hit.y << ", z=" << hit.z
                  << ", new charge = " << hit.charge << ", previous charge = " << origCharge
                  << " smeared threshold = " << actualThreshold << "(delta = " << smear << ")" << endmsg;
    }
}

/**
 * Digitizes the charge.
 * Discretization based on number of bits and bin width scheme.
 */
void MuonCVXDDigitiser::ChargeDigitizer(PixelHitMap& pixels, const HitState& state) const
{
    verbose() << "Charge discretization" << endmsg;

    float minThreshold = m_threshold;
    float maxThreshold = m_chargeMax;
    //int split = 0.3; -- future use
    int numBins = pow(2, m_ChargeDigitizeNumBits.value()) - 1;
    double discCharge = -999;
    for (auto& [pixelID, hit] : pixels) {
        float origCharge = hit.charge;
        discCharge = origCharge;

        switch (m_ChargeDigitizeBinning) {
            case 0: { // uniform binning
                if (origCharge < 1.0) break;
                float binWidth = (maxThreshold - minThreshold) / (numBins);
                if (origCharge < binWidth) discCharge = (minThreshold + binWidth) / 2;
                else if (origCharge > maxThreshold) discCharge = (maxThreshold - binWidth / 2);
                else discCharge = ((ceil((origCharge - binWidth) / binWidth) * binWidth) * 2 + binWidth) / 2;
                break;
            }
            case 1: { // variable binning
                if (origCharge < 1.0) break;
                const std::vector<double>& bins = m_DigitizedBins[state.currentLayer];
                int binVal = -1;
                for (unsigned int idx = 0; idx < bins.size() - 1; idx++) {
                    if (bins[idx + 1] > origCharge) {
                        binVal = idx;
                        break;
                    }
                }
                if (binVal < 0) discCharge = (bins[bins.size() - 2] + bins[bins.size() - 1]) / 2;
                else discCharge = (bins[binVal] + bins[binVal + 1]) / 2;
                break;
            }
        }
        hit.charge = discCharge;
        verbose() << pixelID << ": x=" << hit.x << ", y=" << hit.y << ", z=" << hit.z
                  << ", new charge = " << hit.charge << ", previous charge = " << origCharge
                  << ", number of bits = " << m_ChargeDigitizeNumBits.value()
                  << ", binning scheme = " << m_ChargeDigitizeBinning.value() << endmsg;
    }
}

/**
 * Apply effective measurement resolution.
 * TODO: Right now assuming completely uncorrelated resolution across pixels, will need to divide into:
 * - correlated across pixels, uncorrelated across clusters
 * - correlated within the event, un-correlate
*/
void MuonCVXDDigitiser::TimeSmearer(PixelHitMap& pixels, const HitState& state, EventRandom& random) const
{
    verbose() << "Adding resolution effect to timing measurements" << endmsg;

    // The resolution depends on the layer, not on the individual hit, so it is
    // evaluated once here rather than inside the loop below.
    double sigma_total = m_timeSmearingSigma;

    if (m_timeSmearingModel == 2)
    {
        const double thickness = m_geo.layerThickness[state.currentLayer];
        // -- Realistic timing in planar sensors application default values: -- //
        double t_riseDefault = (8.8 * thickness * 1e3 + 152.1) * 1e-3; //[ns]
        double t_rise = (m_t_riseOverride >= 0.) ? m_t_riseOverride.value() : t_riseDefault;
        double sigma_landauDefault   = 0.03 * thickness / 0.05;  // [ns]from sensor thickness & charge deposition fluctuations - 30ps/50microns
        double sigma_timewalkDefault = 0.1 * t_rise; //[ns] t_rise * 0.1
        double sigma_jitterDefault   = (m_electronicNoise * t_rise) / (80000 * thickness);  // [ns] Q_noise/slope in charge over time 80e/micron = 80000e/mm
        double sigma_TDCDefault      = 0.025 / std::sqrt(12);  // [ns] time to digital converter, 25 ps LSB
        double sigma_clockDefault    = 0.005;  // [ns] fixed by clock quality, 5 ps

        double sigma_landau = (m_sigma_landauOverride >= 0.) ? m_sigma_landauOverride.value() : sigma_landauDefault;
        double sigma_timewalk = (m_sigma_timewalkOverride >= 0.) ? m_sigma_timewalkOverride.value() : sigma_timewalkDefault;
        double sigma_jitter = (m_sigma_jitterOverride >= 0.) ? m_sigma_jitterOverride.value() : sigma_jitterDefault;
        double sigma_TDC = (m_sigma_TDCOverride >= 0.) ? m_sigma_TDCOverride.value() : sigma_TDCDefault;
        double sigma_clock = (m_sigma_clockOverride >= 0.) ? m_sigma_clockOverride.value() : sigma_clockDefault;

        sigma_total = std::sqrt(sigma_landau * sigma_landau + sigma_timewalk * sigma_timewalk + sigma_jitter * sigma_jitter
                                + sigma_TDC * sigma_TDC + sigma_clock * sigma_clock);
    }

    verbose() << "sigma_total: " << sigma_total << endmsg;
    if (sigma_total <= 0.) return;

    for (auto& [pixelID, hit] : pixels)
    {
        double delta = random.gauss(0., sigma_total);

        hit.time = hit.time + delta;
        verbose() << pixelID << ": x=" << hit.x << ", y=" << hit.y << ", z=" << hit.z
                  << ", time = " << hit.time << "(delta = " << delta << ")" << endmsg;
    }
}

/**
 * Digitizes the time information.
 * Discretization based on number of bits and bin width scheme.
 */
void MuonCVXDDigitiser::TimeDigitizer(PixelHitMap& pixels) const
{
    verbose() << "Time discretization" << endmsg;

    if (m_TimeDigitizeBinning != 0) {
        ++m_invalidTimeBinning;
        return;
    }

    const int numBins = pow(2, m_TimeDigitizeNumBits.value()) - 1;
    // uniform binning
    const float binWidth = m_timeMax / numBins;
    double discTime;
    for (auto& [pixelID, hit] : pixels)
    {
        float origTime = hit.time;
        if (origTime < binWidth) discTime = binWidth / 2;
        else if (origTime > m_timeMax) discTime = m_timeMax - binWidth / 2;
        else discTime = ((ceil((origTime - binWidth) / binWidth) * binWidth) * 2 + binWidth) / 2;
        hit.time = discTime;
        verbose() << pixelID << ": x=" << hit.x << ", y=" << hit.y << ", z=" << hit.z
                  << ", new time = " << hit.time << ", previous time = " << origTime << endmsg;
    } //end loop over pixel cells
}

/**
 * Emulates reconstruction of Tracker Hit
 * Tracker hit position is reconstructed as weighted average of edge pixels.
 * The position is corrected for Lorentz shift.
 * Time is the arithmetic average of constituents.
 * Returns false if no pixel has charge above threshold.
 */
bool MuonCVXDDigitiser::ReconstructTrackerHit(const PixelHitMap& pixels, const HitState& state, ClusterHit& cluster) const
{
    double pos[3] = {0, 0, 0};

    double minX = 99999999;
    double minY = 99999999;
    double maxX = -99999999;
    double maxY = -99999999;

    double charge = 0; // total cluster charge
    unsigned int size = 0; // number of pixels in the cluster
    unsigned int edge_size_minx = 0; //number of pixels at the lower edge of cluster in x direction
    unsigned int edge_size_miny = 0; //number of pixels at the lower edge of cluster in y direction
    unsigned int edge_size_maxx = 0; //number of pixels at the upper edge of cluster in x direction
    unsigned int edge_size_maxy = 0; //number of pixels at the upper edge of cluster in y direction

    verbose() << "Creating reconstructed cluster" << endmsg;
    double time = 0; //average time

    /* Get extreme positions, currently only implemented for barrel */
    /* Calculate the mean */
    for (const auto& [pixelID, hit] : pixels)
    {
        //check for non-zero value (pixels below threshold have already been set to zero)
        if (hit.charge < 1.0) continue;

        size += 1;
        time += hit.time;
        charge += hit.charge;
        verbose() << pixelID << ": Averaging position, x=" << hit.x << ", y=" << hit.y << ", weight(EDep)=" << hit.charge << endmsg;

        // calculate min x, min y, max x, max y
        if (hit.x < minX) minX = hit.x;
        if (hit.y < minY) minY = hit.y;
        if (hit.x > maxX) maxX = hit.x;
        if (hit.y > maxY) maxY = hit.y;
    }

    if (not (charge > 0.)) return false;

    // Loop over all pixel hits again, find the pixels on the 4 extreme edges
    for (const auto& [pixelID, hit] : pixels) {
        if (hit.charge < 1.0) continue; // ignore pixels below threshold

        if (hit.x == minX) edge_size_minx += 1;
        if (hit.y == minY) edge_size_miny += 1;
        if (hit.x == maxX) edge_size_maxx += 1;
        if (hit.y == maxY) edge_size_maxy += 1;
    }

    /* Calculate mean x and y by weighted ave:
    x_reco = ( (x_max * max edge size) + (x_min * min edge size) ) / (min edge size + max edge size) */
    pos[0] = ((minX * edge_size_minx) + (maxX * edge_size_maxx)) / (edge_size_minx + edge_size_maxx);
    pos[1] = ((minY * edge_size_miny) + (maxY * edge_size_maxy)) / (edge_size_miny + edge_size_maxy);

    cluster.eDep = (charge / m_electronsPerKeV) * dd4hep::keV;

    const double halfThickness = m_geo.layerHalfThickness[state.currentLayer];
    verbose() << "Edge sizes, minx, maxx, miny, maxy: " << edge_size_minx << ", " << edge_size_maxx << ", "
              << edge_size_miny << ", " << edge_size_maxy << "\n"
              << "Position: x = " << pos[0] << " + " << halfThickness * m_tanLorentzAngleX << "(LA-correction)";
    pos[0] -= halfThickness * m_tanLorentzAngleX;
    verbose() << " = " << pos[0]
              << "; y = " << pos[1] << " + " << halfThickness * m_tanLorentzAngleY << "(LA-correction)";
    pos[1] -= halfThickness * m_tanLorentzAngleY;
    verbose() << " = " << pos[1];

    for (int i = 0; i < 3; ++i) cluster.position[i] = pos[i];
    time /= size;
    cluster.time = time;
    verbose() << ", time (ns) = " << time << endmsg;

    return true;
}

/** Function transforms local coordinates in the ladder
 * into global coordinates
 */
void MuonCVXDDigitiser::TransformToLab(const std::uint64_t cellID, const double* xLoc, double* xLab) const
{
    // Use SurfaceManager to calculate global coordinates
    verbose() << "Cell ID of Hit (used for transforming to lab coords)" << cellID << endmsg;
    const ISurface* surf = m_geo.surfaceMap->find(cellID)->second;
    Vector2D oldPos(xLoc[0] * dd4hep::mm, xLoc[1] * dd4hep::mm);
    Vector3D lv = surf->localToGlobal(oldPos);
    // Store local position in mm
    for (int i = 0; i < 3; i++)
        xLab[i] = lv[i] / dd4hep::mm;
}

/**
 * Function calculates position in pixel matrix based on the
 * local coordinates of point in the sensor.
 */
void MuonCVXDDigitiser::TransformXYToCellID(double x, double y, int& ix, int& iy, const HitState& state) const
{
    // Shift all of L/2 so that all numbers are positive
    iy = int((y + state.sensorHalfLengthV) / m_pixelSizeY);
    ix = int((x + state.sensorHalfLengthU) / m_pixelSizeX);
}

/**
 Function calculates position in the local frame
 based on the index of pixel in the sensor.
*/
void MuonCVXDDigitiser::TransformCellIDToXY(int ix, int iy, double& x, double& y, const HitState& state) const
{
    // Put the point in the cell center
    y = ((0.5 + double(iy)) * m_pixelSizeY) - state.sensorHalfLengthV;
    x = ((0.5 + double(ix)) * m_pixelSizeX) - state.sensorHalfLengthU;
}

// The tolerance keeps a sensor length that is a multiple of the pixel size from getting an extra pixel from rounding
int MuonCVXDDigitiser::GetPixelsInaColumn(const HitState& state) const //SP: why columns!?! I would have guess row..
{
    return std::ceil(2 * state.sensorHalfLengthU / m_pixelSizeX - 1e-6);
}

int MuonCVXDDigitiser::GetPixelsInaRow(const HitState& state) const
{
    return std::ceil(2 * state.sensorHalfLengthV / m_pixelSizeY - 1e-6);
}

void MuonCVXDDigitiser::PrintGeometryInfo() const
{
    std::ostringstream out;
    out << "Pixel size X: " << m_pixelSizeX.value() << "\n"
        << "Pixel size Y: " << m_pixelSizeY.value() << "\n"
        << "Electrons per KeV: " << m_electronsPerKeV.value() << "\n";
    m_geo.print(out);
    info() << out.str() << endmsg;
}

//=============================================================================
// Sample charge from 1 / n^2 distribution.
//=============================================================================
double MuonCVXDDigitiser::randomTail(const double qmin, const double qmax, EventRandom& random) const
{
    const double offset = 1. / qmax;
    const double range  = (1. / qmin) - offset;
    const double u      = offset + random.flat() * range;
    return 1. / u;
}

int MuonCVXDDigitiser::layerMapping(int id) const
{
    int mappedLayerID = m_geo.layerIndex(id);
    if (mappedLayerID < 0 || mappedLayerID >= m_geo.numberOfLayers) {
        ++m_unmappedLayer;
        debug() << id << " not found in the LayerIDs vector." << endmsg;
        return -1;
    }
    verbose() << "Mapped ID of layer " << id << " is: " << mappedLayerID << endmsg;
    return mappedLayerID;
}
