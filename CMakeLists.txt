########################################################
# cmake file for building MuonCVXDDigitiser Package
########################################################
cmake_minimum_required(VERSION 3.23)

project(MuonCVXDDigitiser LANGUAGES CXX)

# please keep this layout for version setting, used by the automatic tagging script
set(${PROJECT_NAME}_VERSION_MAJOR 0)
set(${PROJECT_NAME}_VERSION_MINOR 3)
set(${PROJECT_NAME}_VERSION_PATCH 0)

set(${PROJECT_NAME}_VERSION "${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}.${${PROJECT_NAME}_VERSION_PATCH}")

### DEPENDENCIES ############################################################

find_package(podio REQUIRED)
find_package(EDM4HEP REQUIRED)
find_package(DD4hep REQUIRED COMPONENTS DDRec)
find_package(Gaudi REQUIRED)
find_package(k4FWCore REQUIRED)
find_package(ROOT REQUIRED COMPONENTS Core Hist)
find_package(CLHEP REQUIRED)
find_package(GSL REQUIRED)

include(cmake/Key4hepConfig.cmake)

include(GNUInstallDirs)

### LIBRARIES ###############################################################

add_subdirectory(MuonCVXDDigitiser)

### TESTS ###################################################################

include(CTest)
if(BUILD_TESTING)
  add_subdirectory(test)
endif()

### DOCUMENTATION ###########################################################

option(INSTALL_DOC "Set to OFF to skip build/install Documentation" OFF)

if(INSTALL_DOC)
  add_subdirectory(doc)
endif()

### PACKAGE CONFIG ##########################################################

include(cmake/CreateProjectConfig.cmake)
