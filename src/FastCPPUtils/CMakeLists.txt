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

set(SOURCES
	FastCPPUtils.h rpoly.c math_wrappers.h math_wrappers.cpp
	)

# Library
add_library(FastCPPUtils STATIC ${SOURCES})
set_target_properties(FastCPPUtils PROPERTIES LINKER_LANGUAGE CXX)
if(UNIX)
  target_compile_options(FastCPPUtils PRIVATE -fPIC)
endif()

