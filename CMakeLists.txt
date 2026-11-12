cmake_minimum_required(VERSION 3.10)
project(geodetic_network_adjustment CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(gna src/main.cpp src/GeoNetwork.cpp)
target_include_directories(gna PRIVATE include)

# Eigen (header-only). Uses an installed Eigen3 package if found,
# otherwise pass -DEIGEN3_INCLUDE_DIR=/path/to/eigen
find_package(Eigen3 3.3 QUIET NO_MODULE)
if(TARGET Eigen3::Eigen)
  target_link_libraries(gna PRIVATE Eigen3::Eigen)
else()
  find_path(EIGEN3_INCLUDE_DIR Eigen/Dense PATH_SUFFIXES eigen3)
  if(NOT EIGEN3_INCLUDE_DIR)
    message(FATAL_ERROR "Eigen not found. Install Eigen 3 or set -DEIGEN3_INCLUDE_DIR=/path/to/eigen")
  endif()
  target_include_directories(gna PRIVATE ${EIGEN3_INCLUDE_DIR})
endif()
