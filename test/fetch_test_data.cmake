# Fetch the test samples into DEST: copied from LOCAL_DIR if it has them, otherwise downloaded
# from URL. Every file is checked against its SHA-256.
#
# cmake -DDEST=<dir> -DURL=<release download URL> [-DLOCAL_DIR=<dir>] -P fetch_test_data.cmake

set(samples
  "mu_gun_detailed.slcio=48ce8f8c6f65058866a74ac7f4ddc12cf9fb3487bf4bd84fa6eb5eb15f503d31"
  "mu_gun_nodetail.slcio=126f992d2d4fe8f3a018d78bb7f0decd72525d43cebebd294438df9bf9ff161e"
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
