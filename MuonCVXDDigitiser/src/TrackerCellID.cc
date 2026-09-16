#include "TrackerCellID.h"

TrackerCellID::TrackerCellID(const std::string& encoding) :
    m_encoding(encoding),
    m_coder(encoding),
    m_system(m_coder.index("system")),
    m_side(m_coder.index("side")),
    m_layer(m_coder.index("layer")),
    m_module(m_coder.index("module")),
    m_sensor(m_coder.index("sensor"))
{}

TrackerCellID::CellID TrackerCellID::encode(FieldID system, FieldID side, FieldID layer,
                                            FieldID module, FieldID sensor) const
{
    CellID cellID = 0;
    m_coder.set(cellID, m_system, system);
    m_coder.set(cellID, m_side, side);
    m_coder.set(cellID, m_layer, layer);
    m_coder.set(cellID, m_module, module);
    m_coder.set(cellID, m_sensor, sensor);
    return cellID;
}

TrackerCellID::CellID TrackerCellID::withSensor(CellID cellID, FieldID sensor) const
{
    m_coder.set(cellID, m_sensor, sensor);
    return cellID;
}
