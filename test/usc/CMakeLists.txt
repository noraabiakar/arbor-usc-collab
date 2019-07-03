# Build mechanisms used solely in unit tests.

set(usc_mechanisms
    borgka
    cagk
    cat
    ccanl
    ccanlrev
    gskch
    ichan2
    lca
    nca
    constant_ionic_current
)

include(${PROJECT_SOURCE_DIR}/mechanisms/BuildModules.cmake)

set(external_modcc)
if(ARB_WITH_EXTERNAL_MODCC)
    set(external_modcc MODCC ${modcc})
endif()
set(usc_mech_dir ${CMAKE_CURRENT_BINARY_DIR}/mechanisms)

build_modules(
    ${usc_mechanisms}
    SOURCE_DIR mod
    DEST_DIR "${usc_mech_dir}"
    ${external_modcc}
    MODCC_FLAGS -t cpu -t gpu ${ARB_MODCC_FLAGS} -N usc
    GENERATES .hpp _cpu.cpp _gpu.cpp _gpu.cu
    TARGET build_usc_mods
)

set(usc_mech_sources)
foreach(mech ${usc_mechanisms})
    list(APPEND usc_mech_sources ${usc_mech_dir}/${mech}_cpu.cpp)
    if(ARB_WITH_CUDA)
        list(APPEND usc_mech_sources ${usc_mech_dir}/${mech}_gpu.cpp)
        list(APPEND usc_mech_sources ${usc_mech_dir}/${mech}_gpu.cu)
    endif()
endforeach()

# Test sources

set(usc_sources
    usc_catalogue.cpp
    usc_test.cpp
    test_ccanl.cpp
)

add_executable(usc-test EXCLUDE_FROM_ALL ${usc_sources} ${usc_mech_sources})
add_dependencies(usc-test build_usc_mods)

target_compile_options(usc-test PRIVATE ${ARB_CXXOPT_ARCH})
target_include_directories(usc-test PRIVATE "${CMAKE_CURRENT_BINARY_DIR}")
target_link_libraries(usc-test PRIVATE gtest arbor arbor-private-headers arbor-sup)
