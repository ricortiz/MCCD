
project(MCCD)
cmake_minimum_required(VERSION 2.8.10)

# FindOpenMP
find_package(OpenMP)
if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

option(BUILD_SAMPLE_APP "Build sample application" OFF)

include_directories(Geometry)
include_directories(Srcs)
add_subdirectory(Srcs)

if(BUILD_SAMPLE_APP)
  add_subdirectory(Sample)
endif(BUILD_SAMPLE_APP)
