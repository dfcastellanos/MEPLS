cmake_minimum_required(VERSION 2.8.8)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall")

# create executable
set(TARGET run_sim)
project(${TARGET})
add_executable(${TARGET} "step_2.cc")

# Include MEPLS
set(MEPLS_DIR $ENV{MEPLS_DIR})
set(MEPLS_TUTORIAL_DIR "${MEPLS_DIR}/tutorial")

set(MEPLS_INCLUDE_DIR "${MEPLS_DIR}/inc")
set(MEPLS_TUTORIAL_INCLUDE_DIR "${MEPLS_TUTORIAL_DIR}/inc")

include_directories(${MEPLS_INCLUDE_DIR} ${MEPLS_TUTORIAL_INCLUDE_DIR})
target_include_directories(${TARGET} PUBLIC ${MEPLS_INCLUDE_DIR} ${MEPLS_TUTORIAL_INCLUDE_DIR})

# Link to deal.II
find_package(deal.II 8.0 QUIET HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR})

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

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")


