project(csm)
cmake_minimum_required(VERSION 2.8)
set (BASE_DIR ${PROJECT_SOURCE_DIR}/../)

# Default build type is Release (taken from here: http://www.cmake.org/pipermail/cmake/2008-March/020347.html)
if( NOT CMAKE_BUILD_TYPE )
  set( CMAKE_BUILD_TYPE Release CACHE STRING
       "Choose the type of build, options are: None Debug Release RelWithDebInfo
MinSizeRel."
       FORCE )
  message("No build type requested, defaulting to Release")
endif()


# Open Babel
set (OpenBabel_DIR ${BASE_DIR}/openbabel-files)
if(UNIX)
    set (OpenBabel_PLATFORM_DIR ${OpenBabel_DIR}/unix)
    include_directories(${OpenBabel_PLATFORM_DIR}/include/openbabel-2.0)
    find_library(OpenBabel_LIBS NAMES openbabel HINTS ${OpenBabel_PLATFORM_DIR}/lib)
else()
    set (OpenBabel_PLATFORM_DIR ${OpenBabel_DIR}/Windows)
    include_directories(${OpenBabel_PLATFORM_DIR}/include/)
    find_library(OpenBabel_LIBS_Debug NAMES openbabel-2 HINTS ${OpenBabel_PLATFORM_DIR}/lib/x86/Debug)
    find_library(OpenBabel_LIBS_Release NAMES openbabel-2 HINTS ${OpenBabel_PLATFORM_DIR}/lib/x86/Release)
    set(OpenBabel_LIBS debug ${OpenBabel_LIBS_Debug} optimized ${OpenBabel_LIBS_Release})
    message(${OpenBabel_LIBS})
endif()
     
include_directories(${OpenBabel_DIR}/include)

# Eigen
include_directories(${BASE_DIR}/include)

# Boost
find_package(Boost 1.54 REQUIRED log log_setup system thread filesystem date_time)
include_directories(${Boost_INCLUDE_DIRS})
add_definitions(-DBOOST_ALL_DYN_LINK)
add_definitions(-DBOOST_CHRONO_HEADER_ONLY) # See here http://stackoverflow.com/a/13824645/871910

# C++11
if(UNIX)
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif()

# CSM Specific flags
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_ORDER2_PERMUTER")
if(UNIX)
	set (EXTRA_LIBS ${EXTRA_LIBS} pthread)
endif()

# Executable
add_executable(csm Molecule.cpp parseFunctions.cpp permuter.cpp groupPermuter.cpp elements.cpp 
	logging.cpp mainRot.cpp rpoly.c babelAdapter.cpp math_wrappers.cpp drand48.cpp math_utils.cpp calculations.cpp 
	options.cpp PrintOuts.cpp)
target_link_libraries(csm ${OpenBabel_LIBS} ${Boost_LIBRARIES} ${EXTRA_LIBS})
