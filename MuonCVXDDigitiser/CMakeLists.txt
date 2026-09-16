########################################################
# Physics and helper code, independent of the Gaudi algorithms
########################################################

# Sources are added here as they are ported away from Marlin/LCIO
set(core_sources
  src/G4UniversalFluctuation.cc
  src/LayerGeometry.cc
  src/MyG4UniversalFluctuationForSi.cc
  src/TrackerCellID.cc
)

add_library(MuonCVXDDigitiserCore SHARED ${core_sources})

target_include_directories(MuonCVXDDigitiserCore PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(MuonCVXDDigitiserCore PUBLIC
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
  LINK
    MuonCVXDDigitiserCore
    Gaudi::GaudiKernel
    k4FWCore::k4FWCore
    k4FWCore::k4Interface
    EDM4HEP::edm4hep
    DD4hep::DDRec
    GSL::gsl
)
