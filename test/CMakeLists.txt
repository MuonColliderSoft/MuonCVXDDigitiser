########################################################
# Tests
########################################################

add_executable(test_shared_infrastructure test_shared_infrastructure.cpp)
target_link_libraries(test_shared_infrastructure PRIVATE MuonCVXDDigitiserCore)

if(DEFINED ENV{K4GEO})
  add_test(NAME test_shared_infrastructure
    COMMAND test_shared_infrastructure $ENV{K4GEO}/MuColl/MAIA/compact/MAIA_v0/MAIA_v0.xml)
else()
  message(WARNING "K4GEO is not set: geometry-based tests are disabled")
endif()
