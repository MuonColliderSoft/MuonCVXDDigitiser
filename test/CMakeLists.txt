########################################################
# Unit tests
########################################################

# Choice of the charge digitization bins, which needs no geometry
add_executable(test_charge_bins test_charge_bins.cpp)
target_link_libraries(test_charge_bins PRIVATE MuonCVXDDigitiserCore)
add_test(NAME test_charge_bins COMMAND test_charge_bins)

add_executable(test_sensors test_sensors.cpp)
target_link_libraries(test_sensors PRIVATE MuonCVXDDigitiserCore)
add_test(NAME test_sensors COMMAND test_sensors)

add_executable(test_shared_infrastructure test_shared_infrastructure.cpp)
target_link_libraries(test_shared_infrastructure PRIVATE MuonCVXDDigitiserCore)

if(NOT DEFINED ENV{K4GEO})
  message(WARNING "K4GEO is not set: the tests that need the MAIA_v0 geometry are disabled")
  return()
endif()

add_test(NAME test_shared_infrastructure
  COMMAND test_shared_infrastructure $ENV{K4GEO}/MuColl/MAIA/compact/MAIA_v0/MAIA_v0.xml)

########################################################
# Regression tests of the digitisers on muon-gun samples
#
# The samples are attached to a release of the repository and fetched by the first test. Set
# MUONCVXD_TEST_DATA_DIR to a folder holding them to work offline.
########################################################

set(MUONCVXD_TEST_DATA_URL
    "https://github.com/MuonColliderSoft/MuonCVXDDigitiser/releases/download/edm4hep-test-inputs"
    CACHE STRING "Where the test samples are downloaded from")
set(MUONCVXD_TEST_DATA_DIR "" CACHE PATH "Local folder with the test samples, used instead of downloading them")

find_package(Python3 REQUIRED COMPONENTS Interpreter)

set(run ${PROJECT_BINARY_DIR}/run)
set(data_dir ${CMAKE_CURRENT_BINARY_DIR}/data)
set(out_dir ${CMAKE_CURRENT_BINARY_DIR}/out)
set(reference_dir ${CMAKE_CURRENT_SOURCE_DIR}/reference)
file(MAKE_DIRECTORY ${out_dir})

add_test(NAME fetch_test_data
  COMMAND ${CMAKE_COMMAND} -DDEST=${data_dir} -DURL=${MUONCVXD_TEST_DATA_URL}
          -DLOCAL_DIR=${MUONCVXD_TEST_DATA_DIR} -P ${CMAKE_CURRENT_SOURCE_DIR}/fetch_test_data.cmake)
set_tests_properties(fetch_test_data PROPERTIES FIXTURES_SETUP TestData)

set(reference_commands "")

# digi_<name>: run the digitiser with a configuration from digi_common.py on a sample.
# Extra arguments are passed to make_options.py (--algorithm, --subdet, --set NAME=VALUE, --threads).
function(add_digi_run name config sample)
  add_test(NAME digi_${name}
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/run_digi.sh ${run}
            ${config} ${data_dir}/${sample}.edm4hep.root ${out_dir}/${name} ${ARGN})
  set_tests_properties(digi_${name} PROPERTIES FIXTURES_REQUIRED TestData FIXTURES_SETUP Digi_${name})
endfunction()

# check_<name>: sanity requirements on the output of digi_<name>, and agreement with
# reference/<name>.json. Extra arguments are passed to check_hits.py (e.g. --layers).
function(add_check_test name)
  add_test(NAME check_${name}
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_hits.py ${out_dir}/${name}.csv
            --reference ${reference_dir}/${name}.json ${ARGN})
  set_tests_properties(check_${name} PROPERTIES FIXTURES_REQUIRED Digi_${name})

  set(reference_commands ${reference_commands}
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_hits.py ${out_dir}/${name}.csv
              --write-reference ${reference_dir}/${name}.json ${ARGN}
      PARENT_SCOPE)
endfunction()

# digi_<name> and check_<name>. Arguments after CHECK_ARGS are passed to check_hits.py.
function(add_digi_test name config sample)
  cmake_parse_arguments(PARSE_ARGV 3 arg "" "" "CHECK_ARGS")
  add_digi_run(${name} ${config} ${sample} ${arg_UNPARSED_ARGUMENTS})
  add_check_test(${name} ${arg_CHECK_ARGS})
  set(reference_commands ${reference_commands} PARENT_SCOPE)
endfunction()

# digi_<name> with MuonCVXDRealDigitiser, on the barrel sub-detectors: it does not support endcaps.
# Its hits are linked to all their SimTrackerHits, so mean_pixels is the number of SimTrackerHits
# per cluster.
function(add_real_digi_run name config sample)
  add_digi_run(${name} ${config} ${sample} --algorithm MuonCVXDRealDigitiser
               --subdet VXDBarrel --subdet ITBarrel --subdet OTBarrel ${ARGN})
endfunction()

function(add_real_digi_test name config sample)
  cmake_parse_arguments(PARSE_ARGV 3 arg "" "" "CHECK_ARGS")
  add_real_digi_run(${name} ${config} ${sample} ${arg_UNPARSED_ARGUMENTS})
  add_check_test(${name} ${arg_CHECK_ARGS})
  set(reference_commands ${reference_commands} PARENT_SCOPE)
endfunction()

# compare_<name>: compare two digitisation outputs (options of compare_hits.py)
function(add_compare_test name a b)
  add_test(NAME compare_${name}
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/compare_hits.py
            ${out_dir}/${a}.csv ${out_dir}/${b}.csv --label ${name} ${ARGN})
  set_tests_properties(compare_${name} PROPERTIES FIXTURES_REQUIRED "Digi_${a};Digi_${b}")
endfunction()

# compare_<name>_edm4hep: every hit, fired pixel and link of two outputs is identical
function(add_identical_edm4hep_test name a b)
  add_test(NAME compare_${name}_edm4hep
    COMMAND ${run} ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/compare_digi.py
            ${out_dir}/${a}.edm4hep.root ${out_dir}/${b}.edm4hep.root)
  set_tests_properties(compare_${name}_edm4hep PROPERTIES FIXTURES_REQUIRED "Digi_${a};Digi_${b}")
endfunction()

# All tracker sub-detectors, default settings and without smearing
add_digi_test(default default mu_gun_detailed)
add_digi_test(nosmear nosmear mu_gun_detailed)

# Consistency of the EDM4hep output: one link per hit to a SimTrackerHit on the same sensor, raw hit
# links to every fired pixel with weights summing to 1, and no pixels when they are not requested
add_test(NAME check_default_edm4hep
  COMMAND ${run} ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_digi.py ${out_dir}/default.edm4hep.root --pixels)
set_tests_properties(check_default_edm4hep PROPERTIES FIXTURES_REQUIRED Digi_default)
add_digi_run(nopixels default mu_gun_detailed --set StoreFiredPixels=False)
add_test(NAME check_nopixels_edm4hep
  COMMAND ${run} ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_digi.py ${out_dir}/nopixels.edm4hep.root)
set_tests_properties(check_nopixels_edm4hep PROPERTIES FIXTURES_REQUIRED Digi_nopixels)

# Hits without momentum and path length (ddsim without --enableDetailedShowerMode)
add_digi_test(nosmear_nodetail nosmear mu_gun_nodetail)

# Same configuration twice, and multithreaded: the random numbers are seeded per event and
# algorithm, so the output must be reproducible whatever the order the events are processed in
add_digi_run(nosmear_repeat nosmear mu_gun_detailed)
add_compare_test(reproducible nosmear nosmear_repeat --expect identical)
add_digi_run(default_mt default mu_gun_detailed --threads 4)
add_identical_edm4hep_test(default_mt default default_mt)

# Charge digitization with the 3-bit table and with uniform binning
add_digi_test(bits3 bits3 mu_gun_detailed)
add_digi_test(uniform uniform mu_gun_detailed)

# A layer missing from LayerIDs: its hits are skipped, the other layers are digitised as usual.
# The last vertex-barrel layer (cell ID layer 6) is left out, so the others keep their geometry.
add_digi_test(skip_layer nosmear mu_gun_detailed --subdet VXDBarrel "--set=VXDBarrel.LayerIDs=0 1 2 4"
              CHECK_ARGS --layers VXDBarrel=0,1,2,4)

# Realistic time smearing, derived from the sensor thickness
add_digi_test(timing_realistic timing_realistic mu_gun_detailed)

# Every parameter of the realistic model is used: with all terms overridden, it gives the same hits
# as the constant model with the same sigma
add_digi_run(timing_overrides timing_overrides mu_gun_detailed)
add_digi_run(timing_sigma100 timing_sigma100 mu_gun_detailed)
add_compare_test(timing_overrides timing_overrides timing_sigma100 --expect identical)

# Two instances of the algorithm in one job with different time digitizations: each puts its
# fired-pixel times on its own bins
add_digi_test(time_bins nosmear mu_gun_detailed --subdet VXDBarrel --subdet VXDEndcap
              --set=VXDBarrel.TimeDigitizeNumBits=4 --set=VXDBarrel.TimeMaximum=10
              --set=VXDEndcap.TimeDigitizeNumBits=6 --set=VXDEndcap.TimeMaximum=5
              CHECK_ARGS --time-bins VXDBarrel=4:10 --time-bins VXDEndcap=6:5)

# Ionisation trail re-simulated from the local direction and EnergyLoss
add_digi_test(resim resim mu_gun_detailed)

# The source of the ionisation trail, on mu_gun_scaled: mu_gun_detailed with the Geant4 deposit
# doubled and the path length halved (make_test_inputs.sh). The re-simulated trail takes neither
# from Geant4, so the hits are identical; the Geant4 trail is anchored to the deposit, so the charge
# doubles (within 10 %: the largest vertex pixel charges saturate the last bin of the 4-bit table,
# ratio 1.86 in the barrel).
add_digi_run(resim_scaled resim mu_gun_scaled)
add_digi_run(nosmear_scaled nosmear mu_gun_scaled)
add_compare_test(resim_ignores_geant4 resim resim_scaled --expect identical)
add_compare_test(geant4_charge_follows_deposit nosmear nosmear_scaled --expect-ratio mean_charge_keV=2:0.10)

# Multiple scattering, with and without hit momentum. For 1-100 GeV muons the scattering in the
# sensor is tiny, so cluster size and charge must stay close to the unscattered ones.
add_digi_test(ms ms mu_gun_detailed)
add_digi_test(ms_nodetail ms mu_gun_nodetail)
add_compare_test(ms_vs_nosmear nosmear ms --max-rel-diff mean_pixels=0.05 --max-rel-diff mean_charge_keV=0.05)
# The sample without hit momentum has 5 events (50-200 hits per sub-detector): the means are only
# known to 3-6 %, hence the looser limits.
add_compare_test(ms_vs_nosmear_nodetail nosmear_nodetail ms_nodetail
                 --max-rel-diff mean_pixels=0.10 --max-rel-diff mean_charge_keV=0.10)

# Multiple scattering with the re-simulated ionisation trail, placed along the kinked trajectory
# for the secondaries: the trail must stay in the sensor, so clusters stay close to the unscattered ones.
# Residuals are only checked against the absolute limits of check_hits.py: their median moves by
# up to ~20% between random sequences on this small sample. The resim test is defined above.
add_digi_test(ms_resim ms_resim mu_gun_detailed)
add_compare_test(ms_resim_vs_resim resim ms_resim --max-rel-diff mean_pixels=0.05 --max-rel-diff mean_charge_keV=0.05)

# MuonCVXDRealDigitiser, with the trivial sensor model (default) and the RD53A front-end model
add_real_digi_test(real default mu_gun_detailed)
add_real_digi_test(real_rd53a default mu_gun_detailed --set SensorType=0)
add_real_digi_run(real_repeat default mu_gun_detailed)
add_compare_test(real_reproducible real real_repeat --expect identical)
add_real_digi_run(real_mt default mu_gun_detailed --threads 4)
add_identical_edm4hep_test(real_mt real real_mt)
# Every layer of the three barrels gives hits, each linked to SimTrackerHits on its own sensor
foreach(barrel VXDBarrel:VertexBarrelCollection ITBarrel:InnerTrackerBarrelCollection
               OTBarrel:OuterTrackerBarrelCollection)
  string(REPLACE ":" ";" barrel ${barrel})
  list(GET barrel 0 prefix)
  list(GET barrel 1 sim_collection)
  add_test(NAME check_real_${prefix}_edm4hep
    COMMAND ${run} ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_real_digi.py
            ${out_dir}/real.edm4hep.root --sim ${sim_collection} --prefix ${prefix})
  set_tests_properties(check_real_${prefix}_edm4hep PROPERTIES FIXTURES_REQUIRED Digi_real)
endforeach()

# The example options files of the package
set(options ${PROJECT_SOURCE_DIR}/MuonCVXDDigitiser/options)
add_test(NAME options_MuonCVXDDigitiser
  COMMAND ${run} k4run ${options}/runMuonCVXDDigitiser.py
          --IOSvc.Input ${data_dir}/mu_gun_detailed.edm4hep.root --IOSvc.Output ${out_dir}/options_digi.edm4hep.root)
add_test(NAME options_MuonCVXDRealDigitiser
  COMMAND ${run} k4run ${options}/runMuonCVXDRealDigitiser.py --stats-file ${out_dir}/options_real_stats.root
          --IOSvc.Input ${data_dir}/mu_gun_detailed.edm4hep.root --IOSvc.Output ${out_dir}/options_real.edm4hep.root)
set_tests_properties(options_MuonCVXDDigitiser PROPERTIES FIXTURES_REQUIRED TestData FIXTURES_SETUP OptionsDigi)
set_tests_properties(options_MuonCVXDRealDigitiser PROPERTIES FIXTURES_REQUIRED TestData FIXTURES_SETUP OptionsReal)
add_test(NAME check_options_MuonCVXDDigitiser
  COMMAND ${run} ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_digi.py ${out_dir}/options_digi.edm4hep.root)
add_test(NAME check_options_MuonCVXDRealDigitiser
  COMMAND ${run} ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_real_digi.py ${out_dir}/options_real.edm4hep.root)
set_tests_properties(check_options_MuonCVXDDigitiser PROPERTIES FIXTURES_REQUIRED OptionsDigi)
set_tests_properties(check_options_MuonCVXDRealDigitiser PROPERTIES FIXTURES_REQUIRED OptionsReal)

# Regenerate every reference from the outputs of the last ctest run:
#   ctest && cmake --build . --target update_test_references
add_custom_target(update_test_references ${reference_commands}
  COMMENT "Writing the test references from the last ctest run")
