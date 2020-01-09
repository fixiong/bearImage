CMAKE_MINIMUM_REQUIRED(VERSION 3.0)
PROJECT(bearimage)

message("${BEAR_DIR}")
find_package(TIFF)
find_package(JPEG)
find_package(Boost COMPONENTS filesystem)

INCLUDE_DIRECTORIES(                          
     include
     ../bear/include
     )

SET(SRC_FILES 
src/utility.cpp
src/possion_solver.cpp
src/possion_stiching.cpp
)

SET(ASSIMILATOR_FILES 
assimilator/main.cpp
)

SET(PREVIEWER_FILES 
previewer/main.cpp
previewer/down_semple.cpp
)

SET(VALIDATOR_FILES 
validator/main.cpp
)

add_definitions(-std=c++14 -O3)

add_library(${PROJECT_NAME} SHARED ${SRC_FILES})

add_executable( assimilator ${ASSIMILATOR_FILES}  ${SRC_FILES})
target_link_libraries( assimilator ${TIFF_LIBRARY} ${Boost_LIBRARIES} )

add_executable( previewer ${PREVIEWER_FILES} )
target_link_libraries( previewer ${TIFF_LIBRARY} ${JPEG_LIBRARY} ${Boost_LIBRARIES} )

add_executable( validator ${(VALIDATOR_FILES} )
target_link_libraries( validator ${TIFF_LIBRARY} )

