#ifndef TrackerCellID_h
#define TrackerCellID_h 1

#include <cstdint>
#include <string>

#include "DDSegmentation/BitFieldCoder.h"

/** Encoder/decoder for tracker cell IDs, replacing the LCIO CellIDEncoder/CellIDDecoder
 *  and LCTrackerCellID helpers.
 *
 *  The encoding string is the one of the geometry (e.g. the DD4hep constant
 *  "GlobalTrackerReadoutID"), which must provide the fields
 *  system, side, layer, module and sensor.
 *  Field indices are resolved once at construction, so decoding is cheap.
 */
class TrackerCellID
{
public:
    using CellID = dd4hep::DDSegmentation::CellID;
    using FieldID = dd4hep::DDSegmentation::FieldID;

    /// Throws std::runtime_error if one of the required fields is missing
    explicit TrackerCellID(const std::string& encoding);

    int system(CellID cellID) const { return m_coder.get(cellID, m_system); }
    int side(CellID cellID)   const { return m_coder.get(cellID, m_side); }
    int layer(CellID cellID)  const { return m_coder.get(cellID, m_layer); }
    int module(CellID cellID) const { return m_coder.get(cellID, m_module); }
    int sensor(CellID cellID) const { return m_coder.get(cellID, m_sensor); }

    CellID encode(FieldID system, FieldID side, FieldID layer, FieldID module, FieldID sensor) const;

    /// Returns a copy of cellID with the sensor field replaced
    CellID withSensor(CellID cellID, FieldID sensor) const;

    const std::string& encodingString() const { return m_encoding; }
    const dd4hep::DDSegmentation::BitFieldCoder& coder() const { return m_coder; }

private:
    std::string m_encoding;
    dd4hep::DDSegmentation::BitFieldCoder m_coder;
    size_t m_system;
    size_t m_side;
    size_t m_layer;
    size_t m_module;
    size_t m_sensor;
};

#endif //TrackerCellID_h
