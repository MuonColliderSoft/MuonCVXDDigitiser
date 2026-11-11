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

#include "ChargeDigitizationBins.h"
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
    // One table per layer, chosen from the sensor thickness of the layer
    m_DigitizedBins.resize(m_geo.numberOfLayers);
    for (int i = 0; i < m_geo.numberOfLayers; ++i)
    {
        m_DigitizedBins[i] = chargeDigitizationBins(m_ChargeDigitizeNumBits, m_geo.type.isVertex, m_geo.layerThickness[i]);
        if (m_DigitizedBins[i].empty())
        {
            warning() << "No charge digitization bins defined for " << m_ChargeDigitizeNumBits.value()
                      << " bits; variable binning is unusable for layer " << i << endmsg;
            continue;
        }
        debug() << "Layer " << i << ": sensor thickness " << m_geo.layerThickness[i]
                << " mm, " << m_DigitizedBins[i].size() << " charge bins, first/last = "
                << m_DigitizedBins[i].front() << "/" << m_DigitizedBins[i].back() << endmsg;
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

    // Tracks running in the sensor plane do not cross it and are left unscattered, as are
    // hits with neither momentum nor MC particle.
    const bool msCrosses = std::fabs(dir[2]) > 1e-12 * std::sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);

    state.msTrajectory.clear();
    state.msTotalScatteringAngle = 0.;
    if (m_doMultipleScattering && q_charge > 0. && state.currentParticleMomentum > 0. && msCrosses) {
        verbose() << "Applying multiple scattering formula" << endmsg;
        //**************************************************************************
        //Multiple scattering implementation based on PDG formula (https://pdg.lbl.gov/2020/reviews/rpp2020-rev-passage-particles-matter.pdf)
        //**************************************************************************
        // Momentum and mass from FindLocalPosition(), which falls back to the MC particle
        // when the hit carries no momentum (p = 0 would make theta_0 infinite)
        double p = state.currentParticleMomentum / dd4hep::GeV; //[GeV/c]
        double mass = state.currentParticleMass / dd4hep::GeV; //[GeV/c^2]
        double beta = p / std::sqrt(p * p + mass * mass);

        double x_0 = 93.7; // [mm] -> radiation length in silicon
        double sensorT = m_geo.layerThickness[layer]; // [mm] -> sensor thickness

        // normalize direction
        double mag = std::sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);
        dir[0] /= mag;
        dir[1] /= mag;
        dir[2] /= mag;

        // Step from the face the particle enters through: z = -halfT when it travels towards
        // +z, z = +halfT when it travels towards -z. The entry point is on the straight line
        // through the hit, the exit point is wherever the scattered track leaves the sensor.
        const double zSign = (dir[2] >= 0.) ? 1. : -1.;
        entry[2] = -zSign * 0.5 * sensorT;
        exit[2]  =  zSign * 0.5 * sensorT;
        for (int i = 0; i < 2; ++i) {
            entry[i] = origPos[i] + dir[i] * (entry[2] - origPos[2]) / dir[2];
        }
        for (int i = 0; i < 3; ++i) pos[i] = entry[i];
        const double dirIn[3] = {dir[0], dir[1], dir[2]};

        // Step through the sensor in a whole number of slices no thicker than
        // MSSliceThickness, so the steps cover exactly the sensor thickness. The small
        // relative tolerance absorbs the float layer thickness (50 um / 5 um evaluates to
        // 10.000000149) while still rounding a genuine 10.4 up to 11.
        double sliceT = (m_msSliceThickness > 0.) ? m_msSliceThickness.value() : sensorT;
        int n_slices = std::max(1, (int)std::ceil(sensorT / sliceT * (1. - 1e-6)));
        double z_segment = sensorT / n_slices;

        // theta_0 is computed once for the full sensor crossing with the PDG formula: summing
        // per-slice values underestimates it by 12-21%, because of the log term
        double pathL_total = sensorT / std::fabs(dir[2]);
        double theta_0_total = (0.0136 / (beta * p)) * q_charge * std::sqrt(pathL_total / x_0)
            * (1 + 0.038 * std::log(pathL_total * std::pow(q_charge, 2) / (x_0 * std::pow(beta, 2)))); //as defined in PDG
        state.msTotalScatteringAngle = theta_0_total;
        // Each slice scatters by theta_0_total / sqrt(n_slices), so the slices add up to theta_0_total
        const double theta_0_slice = theta_0_total / std::sqrt((double)n_slices);

        state.msTrajectory.push_back({pos[0], pos[1], pos[2]}); // entry point
        double cumulativePath = 0.;

        bool crossedSensor = true;
        for (int islice = 0; islice < n_slices; ++islice) {
            double theta_plane_x = random.gauss(0., theta_0_slice);
            double theta_plane_y = random.gauss(0., theta_0_slice);

            const double dirOld[3] = {dir[0], dir[1], dir[2]};

            // Rotate the direction (valid for arbitrary angles, unlike tan/atan):
            // in the x-z plane by theta_plane_x, then in the y-z plane by theta_plane_y
            double dirX_rot = dir[0] * std::cos(theta_plane_x) - dir[2] * std::sin(theta_plane_x);
            double dirZ_rot = dir[0] * std::sin(theta_plane_x) + dir[2] * std::cos(theta_plane_x);
            dir[0] = dirX_rot;
            dir[2] = dirZ_rot;
            dirX_rot = dir[1] * std::cos(theta_plane_y) - dir[2] * std::sin(theta_plane_y);
            dirZ_rot = dir[1] * std::sin(theta_plane_y) + dir[2] * std::cos(theta_plane_y);
            dir[1] = dirX_rot;
            dir[2] = dirZ_rot;

            double mag2 = std::sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);
            dir[0] /= mag2;
            dir[1] /= mag2;
            dir[2] /= mag2;

            // Step along the average of the old and new directions: the scattering is spread
            // over the slice
            const double dirAvg[3] = {0.5 * (dirOld[0] + dir[0]),
                                      0.5 * (dirOld[1] + dir[1]),
                                      0.5 * (dirOld[2] + dir[2])};

            // Each step advances exactly one slice in z, so the trajectory ends on the exit
            // face and stays monotonic in z. Stop if the track turns back or runs along the
            // sensor plane, or once the path exceeds MaxTrackLength.
            double dzAvg = dirAvg[2] * zSign;
            if (dzAvg < 1e-6) { crossedSensor = false; break; }
            double pathL_segment = z_segment / dzAvg;
            if (cumulativePath + pathL_segment > m_maxTrkLen) { crossedSensor = false; break; } // pathL_segment >= true step length

            pos[0] += dirAvg[0] * pathL_segment;
            pos[1] += dirAvg[1] * pathL_segment;
            pos[2] = (islice == n_slices - 1) ? exit[2] : pos[2] + zSign * z_segment;

            state.msTrajectory.push_back({pos[0], pos[1], pos[2]});
            cumulativePath += pathL_segment * std::sqrt(dirAvg[0] * dirAvg[0] + dirAvg[1] * dirAvg[1] + dirAvg[2] * dirAvg[2]);
        }
        if (crossedSensor)
        {
            // Exit point: the end of the trajectory, on the exit face
            for (int i = 0; i < 3; ++i) exit[i] = pos[i];
        }
        else
        {
            // The scattered track turned back or ran along the sensor plane: keep the
            // unscattered straight line rather than a trail that does not cross the sensor.
            debug() << "Scattered track did not cross the sensor, multiple scattering not applied" << endmsg;
            for (int i = 0; i < 3; ++i) dir[i] = dirIn[i];
            for (int i = 0; i < 2; ++i)
                exit[i] = origPos[i] + dir[i] * (exit[2] - origPos[2]) / dir[2];
            state.msTrajectory.clear();
        }
    }
    //end of multiple scattering implementation

    for (int i = 0; i < 3; ++i)
        state.currentLocalPosition[i] = origPos[i];
    debug() << "local position: " << state.currentLocalPosition[0] << ", " << state.currentLocalPosition[1] << ", " << state.currentLocalPosition[2] << endmsg;
    // Local unit direction of travel through the sensor. After the multiple-scattering loop
    // above this is the deflected direction, so the re-simulation chain sees its own kinematics.
    double dirMag = std::sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);
    double u[3] = {dir[0] / dirMag, dir[1] / dirMag, dir[2] / dirMag};
    const double halfT = m_geo.layerHalfThickness[layer];
    const bool crossesPlane = std::fabs(u[2]) > 1e-12;

    // Crossing of the full sensor: the distance from the entry to the exit point, which
    // includes the lateral displacement from multiple scattering. A track running in the
    // sensor plane has no crossing and takes the maximum track length.
    double crossingLength = m_maxTrkLen.value();
    if (crossesPlane)
    {
        const double dx = exit[0] - entry[0];
        const double dy = exit[1] - entry[1];
        const double dz = exit[2] - entry[2];
        crossingLength = std::sqrt(dx * dx + dy * dy + dz * dz);
    }

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
    // measured from the hit position. With multiple scattering the depths z of the points are
    // kept and their x, y are taken from the kinked trajectory below.
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

    debug() << "Number of ionization points: " << state.numberOfSegments
            << ", trail from " << (resimulate ? "re-simulation" : "Geant4")
            << ", G4 EDep = " << hcharge << endmsg;

    // When multiple scattering was applied, place the ionisation points on the kinked trajectory.
    // Skip if the scattering angle is small (< 1 mrad) to avoid unnecessary interpolation
    const std::vector<std::array<double, 3>>& traj = state.msTrajectory;
    const bool useKinkedPath = (m_doMultipleScattering && crossesPlane && traj.size() > 1
                                && state.msTotalScatteringAngle > 1e-3);

    for (int i = 0; i < state.numberOfSegments; ++i)
    {
        double sPath = sLo + (i + 0.5) * geomStep;
        double x = origPos[0] + sPath * u[0];
        double y = origPos[1] + sPath * u[1];
        double z = origPos[2] + sPath * u[2];

        if (useKinkedPath)
        {
            // Same depth as on the straight trail, so the extent set above is kept and the point
            // stays in the sensor; x and y interpolated on the trajectory, monotonic in z.
            const double zDir = (traj.back()[2] >= traj.front()[2]) ? 1. : -1.;
            int k = 0;
            while (k + 2 < (int)traj.size() && (z - traj[k + 1][2]) * zDir > 0.) ++k;
            double dzTraj = traj[k + 1][2] - traj[k][2];
            double t = (std::fabs(dzTraj) > 1e-12) ? (z - traj[k][2]) / dzTraj : 0.5;
            t = std::max(0.0, std::min(1.0, t)); // clamp to the trajectory if it stopped early

            x = traj[k][0] + t * (traj[k + 1][0] - traj[k][0]);
            y = traj[k][1] + t * (traj[k + 1][1] - traj[k][1]);
            z = traj[k][2] + t * dzTraj;
        }
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
