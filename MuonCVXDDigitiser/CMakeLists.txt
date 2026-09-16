########################################################
# Physics and helper code, independent of the Gaudi algorithms
########################################################

set(core_sources
  src/AbstractSensor.cc
  src/DetElemSlidingWindow.cc
  src/FindUnionAlgorithm.cc
  src/G4UniversalFluctuation.cc
  src/HitTemporalIndexes.cc
  src/HKBaseSensor.cc
  src/LayerGeometry.cc
  src/MyG4UniversalFluctuationForSi.cc
  src/PixelDigiMatrix.cc
  src/ShapeProcessingSensor.cc
  src/TrackerCellID.cc
  src/TrivialSensor.cc
)

add_library(MuonCVXDDigitiserCore SHARED ${core_sources})

target_include_directories(MuonCVXDDigitiserCore PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(MuonCVXDDigitiserCore PUBLIC
  Gaudi::GaudiKernel
  EDM4HEP::edm4hep
  DD4hep::DDCore
  DD4hep::DDRec
  CLHEP::CLHEP
  GSL::gsl
)

install(TARGETS MuonCVXDDigitiserCore
  EXPORT ${PROJECT_NAME}Targets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
)

########################################################
# Gaudi algorithms
########################################################

gaudi_add_module(MuonCVXDDigitiserPlugins
  SOURCES
    components/MuonCVXDDigitiser.cpp
    components/MuonCVXDRealDigitiser.cpp
  LINK
    MuonCVXDDigitiserCore
    Gaudi::GaudiKernel
    k4FWCore::k4FWCore
    k4FWCore::k4Interface
    EDM4HEP::edm4hep
    DD4hep::DDRec
    GSL::gsl
    ROOT::Hist
)
