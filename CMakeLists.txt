cmake_minimum_required(VERSION 3.18...3.20)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Release or Debug" FORCE)
endif()

project(EMM_TS07D LANGUAGES Fortran)

enable_testing()

set(datestr 2015_075)

set(ts07d_static_url https://rbspgway.jhuapl.edu/sites/default/files/SpaceWeather/TAIL_PAR.zip)
set(ts07d_variable_url https://rbspgway.jhuapl.edu/models/magneticfieldmodeling/ts07d/coeffs_v02/2015/${datestr}.tgz)

set(ts07d_static_dir "${CMAKE_CURRENT_BINARY_DIR}/ts07d_static_dir/")
set(ts07d_static_arc ${ts07d_static_dir}/tail_par.zip)

set(ts07d_variable_dir "${CMAKE_CURRENT_BINARY_DIR}/ts07d_variable_dir/")
set(ts07d_variable_arc ${ts07d_variable_dir}/ts07d_${datestr}.tgz)

if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  string(APPEND CMAKE_Fortran_FLAGS " -std=legacy")
endif()

add_library(ts07d OBJECT ts07d.f)
add_library(geopack OBJECT geopack.f)

configure_file(ts07d_driver.in.f90 ts07d_driver.f90 @ONLY)
add_executable(ts07d_driver ${CMAKE_CURRENT_BINARY_DIR}/ts07d_driver.f90)
target_link_libraries(ts07d_driver PRIVATE ts07d geopack)

add_test(NAME ts07d_driver:${datestr} COMMAND ts07d_driver 5. 6. 7.)

set_tests_properties(ts07d_driver:${datestr} PROPERTIES
 PASS_REGULAR_EXPRESSION "-25.000   -30.157    -9.554"
 TIMEOUT 15)

# download and extract

# avoid anti-leeching
# https//www.useragentstring.com
set(USER_AGENT "Mozilla/5.0 (Windows NT 10.0; WOW64; rv:77.0) Gecko/20100101 Firefox/77.0")

if(NOT EXISTS ${ts07d_static_dir}/tailamebhr1.par)
  if(NOT EXISTS ${ts07d_static_arc})
    file(DOWNLOAD ${ts07d_static_url} ${ts07d_static_arc}
      HTTPHEADER "User-Agent: ${USER_AGENT}"
      INACTIVITY_TIMEOUT 15
      EXPECTED_HASH SHA256=6a51ecfd50b8966ba5592d8cf6d7a99c8f809bc13d9e76fe744d7a93fd0c4c30)
  endif()
  file(ARCHIVE_EXTRACT INPUT ${ts07d_static_arc} DESTINATION ${ts07d_static_dir})
endif()

if(NOT EXISTS ${ts07d_variable_dir}/${datestr}/${datestr}_00_00.par)
  if(NOT EXISTS ${ts07d_variable_arc})
    file(DOWNLOAD ${ts07d_variable_url} ${ts07d_variable_arc}
      HTTPHEADER "User-Agent: ${USER_AGENT}"
      INACTIVITY_TIMEOUT 15)
  endif()
  file(ARCHIVE_EXTRACT INPUT ${ts07d_variable_arc} DESTINATION ${ts07d_variable_dir}/${datestr})
endif()
