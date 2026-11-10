# Fetch the test samples into DEST: copied from LOCAL_DIR if it has them, otherwise downloaded
# from URL. Every file is checked against its SHA-256.
#
# cmake -DDEST=<dir> -DURL=<release download URL> [-DLOCAL_DIR=<dir>] -P fetch_test_data.cmake

# Converted from the LCIO samples of the lcio-test-inputs release
set(samples
  "mu_gun_detailed.edm4hep.root=b21dca6a38fd924bc8d1822f84af515a11f0db020b8b817a9041ae0415559eb6"
  "mu_gun_nodetail.edm4hep.root=82b01a89da85c3991e45985ee116d6f2fd3bd02e4cb54434626cace20c12562c"
  "mu_gun_scaled.edm4hep.root=d24627d5b88709793aaa06ce5ac99f5c2ff94df13088ad2ec3c4e1e14b6f741d"
)

file(MAKE_DIRECTORY ${DEST})
foreach(sample IN LISTS samples)
  string(REPLACE "=" ";" fields "${sample}")
  list(GET fields 0 name)
  list(GET fields 1 sha)
  set(target ${DEST}/${name})

  if(EXISTS ${target})
    file(SHA256 ${target} have)
    if(have STREQUAL sha)
      continue()
    endif()
  endif()

  if(LOCAL_DIR AND EXISTS ${LOCAL_DIR}/${name})
    message(STATUS "Copying ${name} from ${LOCAL_DIR}")
    file(COPY ${LOCAL_DIR}/${name} DESTINATION ${DEST})
  else()
    message(STATUS "Downloading ${URL}/${name}")
    file(DOWNLOAD ${URL}/${name} ${target} STATUS status TLS_VERIFY ON)
    list(GET status 0 code)
    if(NOT code EQUAL 0)
      file(REMOVE ${target})
      message(FATAL_ERROR "Download of ${name} failed: ${status}")
    endif()
  endif()

  file(SHA256 ${target} have)
  if(NOT have STREQUAL sha)
    message(FATAL_ERROR "${name}: SHA-256 ${have} does not match the expected ${sha}")
  endif()
endforeach()
