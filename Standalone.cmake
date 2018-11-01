set(AMPGEN_CXX ${CMAKE_CXX_COMPILER}  CACHE FILEPATH "This should be the path to compiler (use which c++ for macOS)" )

file(GLOB_RECURSE AMPGEN_SRC src/*)
file(GLOB_RECURSE AMPGEN_HDR AmpGen/*)

if(DEFINED ENV{ROOTSYS})
  list(APPEND CMAKE_MODULE_PATH "$ENV{ROOTSYS}/etc/cmake/")
endif()

find_package(ROOT CONFIG REQUIRED COMPONENTS Minuit2 Matrix MathMore MathCore Gpad Tree Graf)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_TEST_OUTPUT_DIRECTORY    "${CMAKE_BINARY_DIR}/bin/test")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

include(CMakeDependentOption)
include(CMakePrintHelpers)

cmake_print_variables(CMAKE_SOURCE_DIR)

# Default build type from the Kitware Blog
set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
    STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(AMPGEN_DEBUG "AmpGen Debug printout")
option(AMPGEN_TRACE "AmpGen Trace printout")

add_library(AmpGen SHARED ${AMPGEN_SRC} ${AMPGEN_HDR})

target_include_directories(AmpGen
  PUBLIC
  "${CMAKE_SOURCE_DIR}"
  SYSTEM PUBLIC
  "${ROOT_INCLUDE_DIRS}")

target_link_libraries(AmpGen
  PUBLIC
  ${ROOT_LIBRARIES}
  ${CMAKE_DL_LIBS})

target_compile_definitions(AmpGen
  PUBLIC
  "AMPGENROOT_CMAKE=\"${CMAKE_BINARY_DIR}/bin\""
  "AMPGEN_CXX=\"${AMPGEN_CXX}\""
  $<$<BOOL:${AMPGEN_DEBUG}>:DEBUGLEVEL=1>
  $<$<BOOL:${AMPGEN_TRACE}>:TRACELEVEL=1>)

target_compile_options(AmpGen
  PUBLIC
  -Wall -Wextra -Wpedantic -g3
  -Wno-unused-parameter
  -Wno-unknown-pragmas
  $<$<CONFIG:Release>:-Ofast>)

find_package(OpenMP)
if(OpenMP_FOUND OR OpenMP_CXX_FOUND)
  if(NOT TARGET OpenMP::OpenMP_CXX)
    add_library(OpenMP::OpenMP_CXX IMPORTED INTERFACE)
    set_property(TARGET OpenMP::OpenMP_CXX
      PROPERTY INTERFACE_COMPILE_OPTIONS ${OpenMP_CXX_FLAGS})
    # Only works if the same flag is passed to the linker; use CMake 3.9+ otherwise (Intel, AppleClang)
    set_property(TARGET OpenMP::OpenMP_CXX
      PROPERTY INTERFACE_LINK_LIBRARIES ${OpenMP_CXX_FLAGS})
    if(CMAKE_VERSION VERSION_LESS 3.4)
      set_property(TARGET OpenMP::OpenMP_CXX APPEND
        PROPERTY INTERFACE_LINK_LIBRARIES -pthread)
    else()
      find_package(Threads REQUIRED)
      set_property(TARGET OpenMP::OpenMP_CXX APPEND
        PROPERTY INTERFACE_LINK_LIBRARIES Threads::Threads)
    endif()
  endif()
  target_link_libraries(AmpGen PUBLIC OpenMP::OpenMP_CXX)
else()
  message(STATUS "OpenMP not found for CXX, you might have forgotten lb-run ROOT bash or CXX=`which g++` in CERN stack")
endif()


# Default to XROOTD only if on CMT system. Can be overridden with -DAMPGEN_XROOTD=ON
if(DEFINED ENV{CMTCONFIG})
  set(AMPGEN_XROOTD_DEFAULT ON)
else()
  set(AMPGEN_XROOTD_DEFAULT OFF)
endif()

cmake_dependent_option(AMPGEN_XROOTD "Turn on XROOTD discovery" ON "AMPGEN_XROOTD_DEFAULT" OFF)

if(AMPGEN_XROOTD)
  find_library(XROOTD_LIB NAMES libXrdCl.so
    HINTS "/cvmfs/lhcb.cern.ch/lib/lcg/releases/LCG_89/xrootd/4.6.0/$ENV{CMTCONFIG}/lib64")
  target_link_libraries(AmpGen PUBLIC ${XROOTD_LIB})
endif()

file(GLOB_RECURSE applications apps/*.cpp )

foreach( file ${applications} )
  get_filename_component( Executable ${file} NAME_WE )
  cmake_print_variables(Executable)
  add_executable(${Executable} ${file})
  target_link_libraries(${Executable} PUBLIC AmpGen)
endforeach()

# file(GLOB_RECURSE options_files options/*.*)
# execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_BINARY_DIR}/bin")
# foreach(file ${options_files})
#   get_filename_component(OptionFile "${file}" NAME)
#   cmake_print_variables(OptionFile)
#   execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink "${file}" "${CMAKE_BINARY_DIR}/bin/${OptionFile}")
# endforeach()

#Setup CMake to run tests
enable_testing()

find_package(Boost COMPONENTS unit_test_framework REQUIRED)
include_directories (${Boost_INCLUDE_DIRS})

file(GLOB TEST_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} test/*.cpp)

foreach(testSrc ${TEST_SRCS})
  get_filename_component(testName ${testSrc} NAME_WE)
  add_executable(${testName} ${testSrc})
  set_target_properties(${testName} 
    PROPERTIES 
    RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_TEST_OUTPUT_DIRECTORY}" 
    RUNTIME_OUTPUT_DIRECTORY_DEBUG   "${CMAKE_TEST_OUTPUT_DIRECTORY}" 
    RUNTIME_OUTPUT_DIRECTORY         "${CMAKE_TEST_OUTPUT_DIRECTORY}" 
    EXECUTABLE_OUTPUT_DIRECTORY      "${CMAKE_TEST_OUTPUT_DIRECTORY}" 
  )
  target_link_libraries(${testName} ${Boost_LIBRARIES} AmpGen)
  message( "Building test: ${testName} in directory =  ${CMAKE_TEST_OUTPUT_DIRECTORY}" )
  add_test(NAME ${testName} WORKING_DIRECTORY ${CMAKE_TEST_OUTPUT_DIRECTORY} COMMAND ${CMAKE_TEST_OUTPUT_DIRECTORY}/${testName} )
endforeach(testSrc)
