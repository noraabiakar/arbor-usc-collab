cmake_minimum_required(VERSION 3.9)
project(arbor-sonata LANGUAGES CXX)

cmake_policy(SET CMP0074 NEW)
set (CMAKE_CXX_STANDARD 14)

# When arbor PR is in, this will become find_package(arbor REQUIRED MPI)
find_package(MPI REQUIRED CXX)
find_package(arbor REQUIRED)
add_executable(
        basket
        basket.cpp)

target_link_libraries(basket PRIVATE arbor::arbor arbor::arborenv ${HDF5_C_LIBRARIES})
target_include_directories(basket PRIVATE ../common/cpp/include ${MPI_CXX_INCLUDE_PATH})

set_target_properties(basket PROPERTIES OUTPUT_NAME arbor-basket)

install(TARGETS basket DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
