CMAKE_MINIMUM_REQUIRED(VERSION 3.0)
PROJECT(test)

find_package(OpenCV REQUIRED)

message(STATUS "OpenCV library status:")
message(STATUS "    version: ${OpenCV_VERSION}")
message(STATUS "    libraries: ${OpenCV_LIBS}")
message(STATUS "    include path: ${OpenCV_INCLUDE_DIRS}")

message("${OPENCV_LIBRARIES}")


INCLUDE_DIRECTORIES(                          
     include
     ${OPENCV_INCLUDE_DIR}
     )

SET(SRC_FILES 
src/image.cpp
src/possion_solver.cpp
src/possion_stiching.cpp
test/main.cpp
)

add_definitions(-std=c++14)

add_executable(${PROJECT_NAME} ${SRC_FILES})

TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${OpenCV_LIBS})