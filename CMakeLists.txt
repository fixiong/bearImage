CMAKE_MINIMUM_REQUIRED(VERSION 3.0)
PROJECT(bearimage)

message("${BEAR_DIR}")


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

add_definitions(-std=c++14)

add_library(${PROJECT_NAME} SHARED ${SRC_FILES})

