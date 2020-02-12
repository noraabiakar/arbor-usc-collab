cmake_minimum_required(VERSION 3.9)
project(arbor-granule LANGUAGES CXX)
ENABLE_LANGUAGE(C)

cmake_policy(SET CMP0074 NEW)
set (CMAKE_CXX_STANDARD 14)

# When arbor PR is in, this will become find_package(arbor REQUIRED MPI)
find_package(MPI REQUIRED CXX)
find_package(arbor REQUIRED)
find_package(HDF5 REQUIRED)
add_executable(
        granule
        granule.cpp
        ../common/cpp/hdf5_lib.cpp)

target_link_libraries(granule PRIVATE arbor::arbor arbor::arborenv ${HDF5_C_LIBRARIES})
target_include_directories(granule PRIVATE ../common/cpp/include ${HDF5_INCLUDE_DIRS} ${MPI_CXX_INCLUDE_PATH})

set_target_properties(granule PROPERTIES OUTPUT_NAME arbor-granule)

install(TARGETS granule DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
