#include "MuonCVXDRealDigitiser.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <sstream>

#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/Detector.h"
#include "DDRec/Surface.h"

#include "DetElemSlidingWindow.h"
#include "EventRandom.h"
#include "HKBaseSensor.h"
#include "HitTemporalIndexes.h"
#include "TrivialSensor.h"

using dd4hep::rec::SurfaceMap;
using dd4hep::rec::ISurface;
using dd4hep::rec::Vector2D;
using dd4hep::rec::Vector3D;

DECLARE_COMPONENT(MuonCVXDRealDigitiser)

MuonCVXDRealDigitiser::MuonCVXDRealDigitiser(const std::string& name, ISvcLocator* svcLoc) :
    MultiTransformer(name, svcLoc,
                     {KeyValues("CollectionName", {"VertexBarrelCollection"}),
                      KeyValues("EventHeader", {"EventHeader"})},
                     {KeyValues("OutputCollectionName", {"VTXTrackerHits"}),
                      KeyValues("RelationColName", {"VTXTrackerHitRelations"})})
{}

StatusCode MuonCVXDRealDigitiser::initialize()
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
        m_geo = loadLayerGeometry(*m_geoSvc->getDetector(), m_subDetName, m_zSegmented, m_layerIDs);
        m_cellID = std::make_unique<TrackerCellID>(m_geoSvc->constantAsString(m_encodingStringVariable));
    } catch (const std::exception& ex) {
        error() << ex.what() << endmsg;
        return StatusCode::FAILURE;
    }

    if (!m_geo.type.isBarrel) {
        error() << "Only barrel sub-detectors are supported, not " << m_subDetName.value() << endmsg;
        return StatusCode::FAILURE;
    }
    if (m_geo.layerIDs.empty()) {
        m_geo.layerIDs.resize(m_geo.numberOfLayers);
        std::iota(m_geo.layerIDs.begin(), m_geo.layerIDs.end(), 0);
    }
    if ((int)m_geo.layerIDs.size() < m_geo.numberOfLayers) {
        error() << "LayerIDs has " << m_geo.layerIDs.size() << " entries but " << m_subDetName.value() << " has "
                << m_geo.numberOfLayers << " layers" << endmsg;
        return StatusCode::FAILURE;
    }
    if (m_window_size <= 0) {
        error() << "WindowSize must be positive" << endmsg;
        return StatusCode::FAILURE;
    }

    debug() << "Z-segmentation " << (m_geo.zSegmented ? "enabled" : "disabled")
            << " (ZSegmented=" << m_zSegmented.value() << ")" << endmsg;

    PrintGeometryInfo();

    if (m_create_stats)
    {
        double max_histox = std::max(m_pixelSizeX, m_pixelSizeY) * 10;
        m_signal_dHisto = std::make_unique<Histogram>(this, "SignalHitDistance", "Signal Hit offset", Gaudi::Accumulators::Axis<double>{1000, 0., max_histox});
        m_bib_dHisto = std::make_unique<Histogram>(this, "BIBHitDistance", "BIB Hit offset", Gaudi::Accumulators::Axis<double>{1000, 0., max_histox});
        m_signal_cSizeHisto = std::make_unique<Histogram>(this, "SignalClusterSize", "Signal Cluster Size", Gaudi::Accumulators::Axis<double>{1000, 0., 50});
        m_signal_xSizeHisto = std::make_unique<Histogram>(this, "SignalClusterSizeinX", "Signal Cluster Size in x", Gaudi::Accumulators::Axis<double>{1000, 0., 20});
        m_signal_ySizeHisto = std::make_unique<Histogram>(this, "SignalClusterSizeinY", "Signal Cluster Size in y", Gaudi::Accumulators::Axis<double>{1000, 0., 20});
        m_signal_zSizeHisto = std::make_unique<Histogram>(this, "SignalClusterSizeinZ", "Signal Cluster Size in z", Gaudi::Accumulators::Axis<double>{1000, 0., 20});
        m_signal_eDepHisto = std::make_unique<Histogram>(this, "SignalClustereDep", "Signal Cluster Energy (MeV)", Gaudi::Accumulators::Axis<double>{1000, 0., 10e-1});
        m_bib_cSizeHisto = std::make_unique<Histogram>(this, "BIBClusterSize", "BIB Cluster Size", Gaudi::Accumulators::Axis<double>{1000, 0., 50});
        m_bib_xSizeHisto = std::make_unique<Histogram>(this, "BIBClusterSizeinX", "BIB Cluster Size in x", Gaudi::Accumulators::Axis<double>{1000, 0., 20});
        m_bib_ySizeHisto = std::make_unique<Histogram>(this, "BIBClusterSizeinY", "BIB Cluster Size in y", Gaudi::Accumulators::Axis<double>{1000, 0., 20});
        m_bib_zSizeHisto = std::make_unique<Histogram>(this, "BIBClusterSizeinZ", "BIB Cluster Size in z", Gaudi::Accumulators::Axis<double>{1000, 0., 20});
        m_bib_eDepHisto = std::make_unique<Histogram>(this, "BIBClustereDep", "BIB Cluster Energy (MeV)", Gaudi::Accumulators::Axis<double>{1000, 0., 10e-1});
    }

    return StatusCode::SUCCESS;
}

std::tuple<edm4hep::TrackerHitPlaneCollection, edm4hep::TrackerHitSimTrackerHitLinkCollection>
MuonCVXDRealDigitiser::operator()(const edm4hep::SimTrackerHitCollection& STHcol,
                                  const edm4hep::EventHeaderCollection& headers) const
{
    edm4hep::TrackerHitPlaneCollection THcol;
    edm4hep::TrackerHitSimTrackerHitLinkCollection relCol;

    m_nSimHits += STHcol.size();
    if (STHcol.empty())
    {
        m_nRecoHits += 0;
        debug() << "Number of produced hits: " << THcol.size() << endmsg;
        return std::make_tuple(std::move(THcol), std::move(relCol));
    }

    const auto seed = m_uidSvc->getUniqueID(headers, name());
    debug() << "Using seed " << seed << " for event " << headers[0].getEventNumber() << " and run "
            << headers[0].getRunNumber() << endmsg;
    EventRandom random(seed);
    MsgStream& log = msgStream();

    std::size_t RELHISTOSIZE { 10 };
    std::vector<std::size_t> relHisto {};
    relHisto.assign(RELHISTOSIZE, 0);

    HitTemporalIndexes t_index { STHcol, *m_cellID };

    for (int layer = 0; layer < m_geo.numberOfLayers; layer++)
    {
        const int layerID = m_geo.layerIDs[layer];
        for (int ladder = 0; ladder < m_geo.laddersInLayer[layer]; ladder++)
        {
            int num_segment_x = 1;
            int num_segment_y = m_geo.zSegmented ? m_geo.sensorsPerLadder[layer] : 1;

            float m_time = t_index.GetMinTime(layerID, ladder);
            if (m_time == HitTemporalIndexes::MAXTIME)
            {
                verbose() << "Undefined min time for layer " << layerID << " ladder " << ladder << endmsg;
                continue;
            }
            //clock time centered at 0
            float nw = floor(fabs(m_time) / m_window_size);
            float start_time = (m_time >= 0) ? nw * m_window_size : -1 * (nw + 1) * m_window_size;

            std::unique_ptr<AbstractSensor> sensor;
            if (m_sensor_type == 1)
            {
                sensor = std::make_unique<TrivialSensor>(layerID, ladder, num_segment_x, num_segment_y,
                                                         m_geo.layerLadderLength[layer], m_geo.layerLadderWidth[layer],
                                                         m_geo.layerThickness[layer], m_pixelSizeX, m_pixelSizeY,
                                                         *m_cellID, m_geo.detectorID, m_threshold,
                                                         start_time, m_window_size, log);
            }
            else
            {
                sensor = std::make_unique<HKBaseSensor>(layerID, ladder, num_segment_x, num_segment_y,
                                                        m_geo.layerLadderLength[layer], m_geo.layerLadderWidth[layer],
                                                        m_geo.layerThickness[layer], m_pixelSizeX, m_pixelSizeY,
                                                        *m_cellID, m_geo.detectorID, m_threshold, m_fe_slope,
                                                        start_time, m_window_size, log);
            }

            if (sensor->GetStatus() != MatrixStatus::ok)
            {
                if (sensor->GetStatus() == MatrixStatus::pixel_number_error)
                    ++m_pixelNumberError;
                else
                    ++m_segmentNumberError;
                debug() << "Geometry error for layer " << layerID << " ladder " << ladder << endmsg;
                continue;
            }

            DetElemSlidingWindow t_window {
                t_index, *sensor,
                m_window_size, start_time,
                m_tanLorentzAngleX, m_tanLorentzAngleY,
                m_cutOnDeltaRays,
                m_diffusionCoefficient,
                m_electronsPerKeV,
                m_segmentLength,
                m_energyLoss,
                3.0,
                m_maxTrkLen,
                m_deltaEne,
                m_geo.surfaceMap,
                m_geo.zSegmented,
                m_fluctuate,
                random,
                log
            };

            while (t_window.active())
            {
                t_window.process();

                SegmentDigiHitList hit_buffer {};
                sensor->buildHits(hit_buffer);

                for (SegmentDigiHit& digiHit : hit_buffer)
                {
                    SurfaceMap::const_iterator sI = m_geo.surfaceMap->find(digiHit.cellID);
                    if (sI == m_geo.surfaceMap->end())
                    {
                        ++m_missingSurface;
                        debug() << "No surface for cell ID " << digiHit.cellID << endmsg;
                        continue;
                    }
                    const ISurface* surf = sI->second;

                    auto recoHit = THcol.create();
                    recoHit.setEDep((digiHit.charge / m_electronsPerKeV) * dd4hep::keV);

                    double loc_pos[3] = {
                        digiHit.x - m_geo.layerHalfThickness[layer] * m_tanLorentzAngleX,
                        digiHit.y - m_geo.layerHalfThickness[layer] * m_tanLorentzAngleY,
                        0
                    };

                    recoHit.setCellID(digiHit.cellID);

                    // See DetElemSlidingWindow::StoreSignalPoints
                    int segment_id = m_cellID->sensor(digiHit.cellID);
                    float s_offset = sensor->GetSensorCols() * sensor->GetPixelSizeY();
                    s_offset *= (float(segment_id) + 0.5);
                    s_offset -= sensor->GetHalfLength();

                    Vector2D oldPos(loc_pos[0] * dd4hep::mm, (loc_pos[1] - s_offset) * dd4hep::mm);
                    Vector3D lv = surf->localToGlobal(oldPos);

                    recoHit.setPosition({lv[0] / dd4hep::mm, lv[1] / dd4hep::mm, lv[2] / dd4hep::mm});

                    recoHit.setTime(digiHit.time);

                    Vector3D u = surf->u();
                    Vector3D v = surf->v();

                    recoHit.setU({float(u.theta()), float(u.phi())});
                    recoHit.setV({float(v.theta()), float(v.phi())});

                    // ALE Does this make sense??? TO CHECK
                    recoHit.setDu(m_pixelSizeX / std::sqrt(12));
                    recoHit.setDv(m_pixelSizeY / std::sqrt(12));

                    bool sig = false;
                    double minx = 999;
                    double maxx = -999;
                    double miny = 999;
                    double maxy = -999;
                    double minz = 999;
                    double maxz = -999;

                    //All the sim-hits are registered for a given reco-hit
                    for (const auto& [index, st_item] : digiHit.sim_hits)
                    {
                        if (m_create_stats)
                        {
                            const auto& stPos = st_item.getPosition();
                            if (!st_item.isOverlay()) sig = true;
                            minx = std::min(minx, stPos.x);
                            maxx = std::max(maxx, stPos.x);
                            miny = std::min(miny, stPos.y);
                            maxy = std::max(maxy, stPos.y);
                            minz = std::min(minz, stPos.z);
                            maxz = std::max(maxz, stPos.z);

                            const auto& recoPos = recoHit.getPosition();
                            double dist = std::sqrt(std::pow(recoPos.x - stPos.x, 2) + std::pow(recoPos.y - stPos.y, 2)
                                                    + std::pow(recoPos.z - stPos.z, 2));
                            ++(*(st_item.isOverlay() ? m_bib_dHisto : m_signal_dHisto))[dist];
                        }
                        auto t_rel = relCol.create();
                        t_rel.setFrom(recoHit);
                        t_rel.setTo(st_item);
                        t_rel.setWeight(1.0);
                    }

                    m_nSimPerReco += digiHit.sim_hits.size();
                    if (digiHit.sim_hits.size() < RELHISTOSIZE)
                    {
                        relHisto[digiHit.sim_hits.size()]++;
                    }

                    if (m_create_stats)
                    {
                        // cluster size histograms
                        const double eDepMeV = recoHit.getEDep() / dd4hep::MeV;
                        if (!sig)
                        {
                            ++(*m_bib_cSizeHisto)[digiHit.size];
                            ++(*m_bib_xSizeHisto)[maxx - minx];
                            ++(*m_bib_ySizeHisto)[maxy - miny];
                            ++(*m_bib_zSizeHisto)[maxz - minz];
                            ++(*m_bib_eDepHisto)[eDepMeV];
                        } else {
                            ++(*m_signal_cSizeHisto)[digiHit.size];
                            ++(*m_signal_xSizeHisto)[maxx - minx];
                            ++(*m_signal_ySizeHisto)[maxy - miny];
                            ++(*m_signal_zSizeHisto)[maxz - minz];
                            ++(*m_signal_eDepHisto)[eDepMeV];
                        }
                    }

                    if (msgLevel(MSG::VERBOSE))
                    {
                        const auto& recoPos = recoHit.getPosition();
                        verbose() << "Reconstructed pixel cluster for "
                                  << sensor->GetLayer() << ":" << sensor->GetLadder()
                                  << ":" << segment_id << "\n"
                                  << "- global position (x,y,z,t) = " << recoPos.x
                                  << ", " << recoPos.y
                                  << ", " << recoPos.z
                                  << ", " << recoHit.getTime() << "\n"
                                  << "- charge = " << recoHit.getEDep() << endmsg;
                    }
                }
            }

            m_nOffSurface += t_window.get_off_surface();
        }
    }

    m_nRecoHits += THcol.size();
    if (msgLevel(MSG::DEBUG))
    {
        debug() << "Number of produced hits: " << THcol.size() << "\n"
                << "Hit relation histogram:\n";
        std::size_t count = 0;
        for (std::size_t k = 0; k < RELHISTOSIZE; k++)
        {
            debug() << k << " " << relHisto[k] << "\n";
            count += relHisto[k];
        }
        debug() << "> " << THcol.size() - count << endmsg;
    }

    return std::make_tuple(std::move(THcol), std::move(relCol));
}

void MuonCVXDRealDigitiser::PrintGeometryInfo() const
{
    std::ostringstream out;
    out << "Pixel size X: " << m_pixelSizeX.value() << "\n"
        << "Pixel size Y: " << m_pixelSizeY.value() << "\n"
        << "Electrons per KeV: " << m_electronsPerKeV.value() << "\n";
    m_geo.print(out);
    info() << out.str() << endmsg;
}
