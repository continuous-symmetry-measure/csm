project(csm)
cmake_minimum_required(VERSION 2.8)
set (BASE_DIR ${PROJECT_SOURCE_DIR}/../)

# Open Babel
set (OPENBABEL_DIR ${BASE_DIR}/openbabel-files/unix)
include_directories(${BASE_DIR}/include)
include_directories(${OPENBABEL_DIR}/include/openbabel-2.0)
link_directories(${OPENBABEL_DIR}/lib)

# Boost
find_package(Boost 1.54 REQUIRED log log_setup system thread filesystem date_time)
include_directories(${Boost_INCLUDE_DIRS})
add_definitions(-DBOOST_LOG_DYN_LINK)

# C++11
SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# CSM Specific flags
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_ORDER2_PERMUTER")
if(UNIX)
	set (EXTRA_LIBS ${EXTRA_LIBS} pthread)
endif()

# Executable
add_executable(csm Molecule.cpp parseFunctions.cpp permuter.cpp groupPermuter.cpp elements.cpp 
	logging.cpp mainRot.cpp rpoly.c babelAdapter.cpp math_wrappers.cpp drand48.cpp math_utils.cpp calculations.cpp 
	options.cpp PrintOuts.cpp)
target_link_libraries(csm openbabel ${Boost_LIBRARIES} ${EXTRA_LIBS})
