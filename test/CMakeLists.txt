########################################################
# Tests
########################################################

add_executable(test_shared_infrastructure test_shared_infrastructure.cpp)
target_link_libraries(test_shared_infrastructure PRIVATE MuonCVXDDigitiserCore)

add_executable(test_sensors test_sensors.cpp)
target_link_libraries(test_sensors PRIVATE MuonCVXDDigitiserCore)
add_test(NAME test_sensors COMMAND test_sensors)

if(DEFINED ENV{K4GEO})
  add_test(NAME test_shared_infrastructure
    COMMAND test_shared_infrastructure $ENV{K4GEO}/MuColl/MAIA/compact/MAIA_v0/MAIA_v0.xml)
else()
  message(WARNING "K4GEO is not set: geometry-based tests are disabled")
endif()

########################################################
# MuonCVXDDigitiser on a small MAIA muon-gun sample
########################################################

set(options ${PROJECT_SOURCE_DIR}/MuonCVXDDigitiser/options)
set(pixels_args "")
foreach(prefix VXDBarrel VXDEndcap ITBarrel ITEndcap OTBarrel OTEndcap)
  list(APPEND pixels_args --${prefix}Digitiser.StoreFiredPixels True)
endforeach()

if(DEFINED ENV{K4GEO})
  add_test(NAME run_ddsim_maia
    COMMAND ddsim --compactFile $ENV{K4GEO}/MuColl/MAIA/compact/MAIA_v0/MAIA_v0.xml
                  --enableGun --gun.particle mu- --gun.energy 10*GeV --gun.distribution uniform
                  --gun.multiplicity 10 --random.seed 42 --numberOfEvents 10
                  --outputFile sim_maia.edm4hep.root)
  set_tests_properties(run_ddsim_maia PROPERTIES FIXTURES_SETUP SimMAIA)

  add_test(NAME run_MuonCVXDDigitiser
    COMMAND ${PROJECT_BINARY_DIR}/run k4run ${options}/runMuonCVXDDigitiser.py
                  --IOSvc.Input sim_maia.edm4hep.root --IOSvc.Output digi_maia.edm4hep.root)
  add_test(NAME run_MuonCVXDDigitiser_pixels
    COMMAND ${PROJECT_BINARY_DIR}/run k4run ${options}/runMuonCVXDDigitiser.py
                  --IOSvc.Input sim_maia.edm4hep.root --IOSvc.Output digi_maia_pixels.edm4hep.root ${pixels_args})
  add_test(NAME run_MuonCVXDDigitiser_pixels_mt
    COMMAND ${PROJECT_BINARY_DIR}/run k4run ${options}/runMuonCVXDDigitiser.py --threads 4
                  --IOSvc.Input sim_maia.edm4hep.root --IOSvc.Output digi_maia_pixels_mt.edm4hep.root ${pixels_args})
  set_tests_properties(run_MuonCVXDDigitiser PROPERTIES FIXTURES_REQUIRED SimMAIA FIXTURES_SETUP DigiMAIA)
  set_tests_properties(run_MuonCVXDDigitiser_pixels PROPERTIES FIXTURES_REQUIRED SimMAIA FIXTURES_SETUP DigiMAIAPixels)
  set_tests_properties(run_MuonCVXDDigitiser_pixels_mt PROPERTIES FIXTURES_REQUIRED SimMAIA FIXTURES_SETUP DigiMAIAPixelsMT)

  add_test(NAME check_MuonCVXDDigitiser
    COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/check_digi.py digi_maia.edm4hep.root)
  add_test(NAME check_MuonCVXDDigitiser_pixels
    COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/check_digi.py digi_maia_pixels.edm4hep.root --pixels)
  add_test(NAME compare_MuonCVXDDigitiser_mt
    COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/compare_digi.py digi_maia_pixels.edm4hep.root digi_maia_pixels_mt.edm4hep.root)
  set_tests_properties(check_MuonCVXDDigitiser PROPERTIES FIXTURES_REQUIRED DigiMAIA)
  set_tests_properties(check_MuonCVXDDigitiser_pixels PROPERTIES FIXTURES_REQUIRED DigiMAIAPixels)
  set_tests_properties(compare_MuonCVXDDigitiser_mt PROPERTIES FIXTURES_REQUIRED "DigiMAIAPixels;DigiMAIAPixelsMT")
endif()

########################################################
# MuonCVXDRealDigitiser on the same sample
########################################################

if(DEFINED ENV{K4GEO})
  foreach(sensor_type 0 1)
    add_test(NAME run_MuonCVXDRealDigitiser_sensor${sensor_type}
      COMMAND ${PROJECT_BINARY_DIR}/run k4run ${options}/runMuonCVXDRealDigitiser.py --sensor-type ${sensor_type}
                    --stats-file real_stats_sensor${sensor_type}.root
                    --IOSvc.Input sim_maia.edm4hep.root --IOSvc.Output digi_real_sensor${sensor_type}.edm4hep.root)
    add_test(NAME run_MuonCVXDRealDigitiser_sensor${sensor_type}_mt
      COMMAND ${PROJECT_BINARY_DIR}/run k4run ${options}/runMuonCVXDRealDigitiser.py --sensor-type ${sensor_type} --threads 4
                    --IOSvc.Input sim_maia.edm4hep.root --IOSvc.Output digi_real_sensor${sensor_type}_mt.edm4hep.root)
    set_tests_properties(run_MuonCVXDRealDigitiser_sensor${sensor_type} PROPERTIES
      FIXTURES_REQUIRED SimMAIA FIXTURES_SETUP DigiRealSensor${sensor_type})
    set_tests_properties(run_MuonCVXDRealDigitiser_sensor${sensor_type}_mt PROPERTIES
      FIXTURES_REQUIRED SimMAIA FIXTURES_SETUP DigiRealSensor${sensor_type}MT)

    add_test(NAME check_MuonCVXDRealDigitiser_sensor${sensor_type}
      COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/check_real_digi.py digi_real_sensor${sensor_type}.edm4hep.root)
    add_test(NAME compare_MuonCVXDRealDigitiser_sensor${sensor_type}_mt
      COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/compare_digi.py
                      digi_real_sensor${sensor_type}.edm4hep.root digi_real_sensor${sensor_type}_mt.edm4hep.root)
    set_tests_properties(check_MuonCVXDRealDigitiser_sensor${sensor_type} PROPERTIES
      FIXTURES_REQUIRED DigiRealSensor${sensor_type})
    set_tests_properties(compare_MuonCVXDRealDigitiser_sensor${sensor_type}_mt PROPERTIES
      FIXTURES_REQUIRED "DigiRealSensor${sensor_type};DigiRealSensor${sensor_type}MT")
  endforeach()
endif()
