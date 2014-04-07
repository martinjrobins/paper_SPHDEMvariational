cmake_minimum_required(VERSION 2.8)
find_package(VTK REQUIRED)
find_package(Boost 1.50.0 REQUIRED)
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-Wno-deprecated -std=c++0x -fopenmp -D_GLIBCXX_PARALLEL")
#set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-Wno-deprecated -std=c++0x")
set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-undefined")

include_directories(${VTK_INCLUDE_DIRS} ${Boost_INCLUDE_DIR} src Aboria/src Aboria/eigen3)

add_executable (single_particle src/single_particle.cpp)
TARGET_LINK_LIBRARIES(single_particle ${VTK_LIBRARIES})
add_executable (multiple_particle src/multiple_particle.cpp)
TARGET_LINK_LIBRARIES(multiple_particle ${VTK_LIBRARIES})
