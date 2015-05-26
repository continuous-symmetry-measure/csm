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

# Eigen
include_directories(${BASE_DIR}/include)

# Boost
if(WIN32)
  set(Boost_USE_STATIC_LIBS ON)
endif()
find_package(Boost 1.54 REQUIRED log log_setup system thread filesystem date_time)
include_directories(${Boost_INCLUDE_DIRS})
if(UNIX)
  add_definitions(-DBOOST_ALL_DYN_LINK)
endif()
add_definitions(-DBOOST_CHRONO_HEADER_ONLY) # See here http://stackoverflow.com/a/13824645/871910

# CSM Specific flags
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_ORDER2_PERMUTER")
if(UNIX)
	set (EXTRA_LIBS ${EXTRA_LIBS} pthread)
endif()

set(SOURCES
	Molecule.cpp parseFunctions.cpp permuter.cpp groupPermuter.cpp elements.cpp 
	logging.cpp mainRot.cpp rpoly.c math_wrappers.cpp drand48.cpp math_utils.cpp calculations.cpp 
	options.cpp
	calculations.h dmatrix.h drand48.h dvector.h elements.h
	groupPermuter.h logging.h math_utils.h math_wrappers.h Molecule.h options.h
	parseFunctions.h permuter.h 
	)

# Library
add_library(csmlib STATIC ${SOURCES} csmlib.cpp csmlib.h)
set_target_properties(csmlib PROPERTIES LINKER_LANGUAGE CXX)
if(UNIX)
  target_compile_options(csmlib PRIVATE -fPIC)
endif()

