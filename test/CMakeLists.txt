########################################################
# Regression tests of MuonCVXDDigitiser in Marlin
#
# The muon-gun samples are attached to a release of the repository and fetched by the first
# test. Set MUONCVXD_TEST_DATA_DIR to a folder holding them to work offline.
########################################################

set(MUONCVXD_TEST_DATA_URL
    "https://github.com/MuonColliderSoft/MuonCVXDDigitiser/releases/download/lcio-test-inputs"
    CACHE STRING "Where the test samples are downloaded from")
set(MUONCVXD_TEST_DATA_DIR "" CACHE PATH "Local folder with the test samples, used instead of downloading them")

# Unit test of the choice of the charge digitization bins, which needs no geometry
add_executable(test_charge_bins test_charge_bins.cpp ${PROJECT_SOURCE_DIR}/src/ChargeDigitizationBins.cc)
target_include_directories(test_charge_bins PRIVATE ${PROJECT_SOURCE_DIR}/include)
add_test(NAME test_charge_bins COMMAND test_charge_bins)

if(NOT DEFINED ENV{K4GEO})
  message(WARNING "K4GEO is not set: the MuonCVXDDigitiser tests need the MAIA_v0 geometry and are disabled")
  return()
endif()

find_package(Python3 REQUIRED COMPONENTS Interpreter)

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
# Extra arguments are passed to make_steering.py (--subdet, --set NAME=VALUE).
function(add_digi_run name config sample)
  add_test(NAME digi_${name}
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/run_digi.sh $<TARGET_FILE:MuonCVXDDigitiser>
            ${config} ${data_dir}/${sample}.slcio ${out_dir}/${name} ${ARGN})
  set_tests_properties(digi_${name} PROPERTIES FIXTURES_REQUIRED TestData FIXTURES_SETUP Digi_${name})
endfunction()

# digi_<name> as above, then check_<name>: sanity requirements, and agreement with
# reference/<name>.json. Arguments after CHECK_ARGS are passed to check_hits.py.
function(add_digi_test name config sample)
  cmake_parse_arguments(PARSE_ARGV 3 arg "" "" "CHECK_ARGS")
  add_digi_run(${name} ${config} ${sample} ${arg_UNPARSED_ARGUMENTS})

  add_test(NAME check_${name}
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_hits.py ${out_dir}/${name}.csv
            --reference ${reference_dir}/${name}.json ${arg_CHECK_ARGS})
  set_tests_properties(check_${name} PROPERTIES FIXTURES_REQUIRED Digi_${name})

  set(reference_commands ${reference_commands}
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_hits.py ${out_dir}/${name}.csv
              --write-reference ${reference_dir}/${name}.json ${arg_CHECK_ARGS}
      PARENT_SCOPE)
endfunction()

# compare_<name>: compare two digitisation outputs (options of compare_hits.py)
function(add_compare_test name a b)
  add_test(NAME compare_${name}
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/compare_hits.py
            ${out_dir}/${a}.csv ${out_dir}/${b}.csv --label ${name} ${ARGN})
  set_tests_properties(compare_${name} PROPERTIES FIXTURES_REQUIRED "Digi_${a};Digi_${b}")
endfunction()

# All tracker sub-detectors, default settings and without smearing
add_digi_test(default default mu_gun_detailed)
add_digi_test(nosmear nosmear mu_gun_detailed)

# Hits without momentum and path length (ddsim without --enableDetailedShowerMode)
add_digi_test(nosmear_nodetail nosmear mu_gun_nodetail)

# Same configuration twice: the output must be reproducible
add_digi_run(nosmear_repeat nosmear mu_gun_detailed)
add_compare_test(reproducible nosmear nosmear_repeat --expect identical)

# Charge digitization with the 3-bit table and with uniform binning
add_digi_test(bits3 bits3 mu_gun_detailed)
add_digi_test(uniform uniform mu_gun_detailed)

# A layer missing from LayerIDs: its hits are skipped, the other layers are digitised as usual.
# The last vertex-barrel layer (cell ID layer 6) is left out, so the others keep their geometry.
add_digi_test(skip_layer nosmear mu_gun_detailed --subdet VXDBarrel "--set=VXDBarrel.LayerIDs=0 1 2 4"
              CHECK_ARGS --layers VXDBarrel=0,1,2,4)

# Regenerate every reference from the outputs of the last ctest run:
#   ctest && cmake --build . --target update_test_references
add_custom_target(update_test_references ${reference_commands}
  COMMENT "Writing the test references from the last ctest run")
