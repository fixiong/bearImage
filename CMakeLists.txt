CMAKE_MINIMUM_REQUIRED(VERSION 3.0)
PROJECT(bearimage)

message("${BEAR_DIR}")
find_package(TIFF)

INCLUDE_DIRECTORIES(                          
     include
     ../bear/include
     )

SET(SRC_FILES 
src/utility.cpp
src/possion_solver.cpp
src/possion_stiching.cpp
src/possion_wrap.cpp
)

SET(ASSIMILATOR_FILES 
assimilator/main.cpp
)

add_definitions(-std=c++14 -O3)

add_library(${PROJECT_NAME} SHARED ${SRC_FILES})

add_executable( assimilator ${ASSIMILATOR_FILES} )
target_link_libraries( assimilator ${TIFF_LIBRARY} ${PROJECT_NAME} )

