cmake_minimum_required(VERSION 2.8.8)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall")

# create executable
set(TARGET run_sim)
project(${TARGET})
add_executable(${TARGET} main.cc)

# include MEPLS
set(MEPLS_DIR $ENV{MEPLS_DIR})

if (NOT DEFINED MEPLS_DIR)
    message(FATAL_ERROR "\n"
            "*** Could not locate MEPLS ***\n\n"
            "You need to set an environment variable MEPLS_DIR with the path to the MEPLS root directory ($ export MEPLS_DIR=\"/path/to/MEPLS\"), or pass a flag -DMEPLS_DIR=/path/to/MEPLS to cmake.\n"
            )
endif ()

set(MEPLS_INCLUDE_DIR "${MEPLS_DIR}/inc")
include_directories(${MEPLS_INCLUDE_DIR})
target_include_directories(${TARGET} PUBLIC ${MEPLS_INCLUDE_DIR})

# Link to deal.II
find_package(deal.II 9.0 REQUIRED QUIET HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR})

if (NOT ${deal.II_FOUND})
    message(FATAL_ERROR "\n"
            "*** Could not locate deal.II. ***\n\n"
            "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
            "or set an environment variable \"DEAL_II_DIR\" that contains this path."
            )
endif ()

include_directories(${DEAL_II_INCLUDE_DIRS})
set(CMAKE_CXX_FLAGS ${DEAL_II_CXX_FLAGS})
set(CMAKE_CXX_FLAGS_RELEASE ${DEAL_II_CXX_FLAGS_RELEASE})
target_link_libraries(${TARGET} ${DEAL_II_LIBRARIES})

# Link to Hdf5
find_package(HDF5 REQUIRED COMPONENTS CXX) # if it doesn't work, define HDF5_DIR with the path to the directory containing lib/ and include/
target_link_libraries(${TARGET} ${HDF5_LIBRARIES})
target_include_directories(${TARGET} PUBLIC ${HDF5_INCLUDE_DIRS})

# Set DEBUG macro with -DDEFINE_DEBUG=on
option(DEFINE_DEBUG "Build the project using debugging code" OFF)
if (DEFINE_DEBUG)
    add_definitions(-DDEBUG)
    set(CMAKE_BUILD_TYPE Debug)
    message("-- Build type is debug")
else()
    message("-- Build type is release")
endif (DEFINE_DEBUG)

# Set OPENMP macro with -DOPENMP=on
option(OPENMP "Enable openMP parallelization" ON)
if (OPENMP)
    message("-- openMP is enabled...")
    add_definitions(-DOPENMP)
else ()
    message("-- OPENMP is disabled...")
endif (OPENMP)

# maximum gcc optimization and enable openMP
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fopenmp")

# create ./results directory
execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory "results")

# add "make run" command
add_custom_target(run COMMAND ${TARGET} DEPENDS ${TARGET} WORKING_DIRECTORY ${CMAKE_PROJECT_DIR})


