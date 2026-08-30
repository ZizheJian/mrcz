cmake_minimum_required(VERSION 3.10)
project(mrcz)

find_package(OpenMP REQUIRED COMPONENTS CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build shared libraries" FORCE)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(tinyxml2_BUILD_TESTING OFF CACHE BOOL "Build tinyxml2 tests" FORCE)
set(BUILD_SZ3_BINARY OFF CACHE BOOL "Build SZ3 command-line and C API tools" FORCE)
set(SZ3_USE_BUNDLED_ZSTD OFF CACHE BOOL "Use SZ3 bundled Zstd" FORCE)
set(SZ_FIND_DEPS ON CACHE BOOL "Use system dependencies for SZ2" FORCE)
set(BUILD_SZ_EXAMPLES OFF CACHE BOOL "Build SZ2 examples" FORCE)
set(BUILD_TESTS OFF CACHE BOOL "Build SZ2 tests" FORCE)
set(BUILD_OPENMP OFF CACHE BOOL "Build SZ2 internal OpenMP" FORCE)

add_subdirectory(../third_party/parsers/c_codes ${CMAKE_BINARY_DIR}/parsers_build)
add_subdirectory(mrcz_core)
add_subdirectory(../third_party/zstd ${CMAKE_BINARY_DIR}/zstd_build)
add_subdirectory(../third_party/zfp ${CMAKE_BINARY_DIR}/zfp_build)
add_subdirectory(../third_party/tinyxml2 ${CMAKE_BINARY_DIR}/tinyxml2_build)
add_subdirectory(../third_party/SZ3 ${CMAKE_BINARY_DIR}/sz3_build EXCLUDE_FROM_ALL)
add_subdirectory(../third_party/SZ2 ${CMAKE_BINARY_DIR}/sz2_build EXCLUDE_FROM_ALL)

add_library(mrcz_qoz_bridge STATIC qoz_impl/qoz_bridge.cpp)
target_include_directories(mrcz_qoz_bridge SYSTEM PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../third_party/QoZ/include)
target_compile_definitions(mrcz_qoz_bridge PRIVATE QoZ_DEBUG_TIMINGS=0)
target_link_libraries(mrcz_qoz_bridge PRIVATE libzstd_static OpenMP::OpenMP_CXX)

set(SPERR_INSTALL_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../install/sperr_install")
set(MFEM_INSTALL_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../install/mfem_install")

function(add_mrcz_executable target source)
    add_executable(${target} ${source})
    target_include_directories(${target} PRIVATE ${MFEM_INSTALL_PATH}/include)
    target_include_directories(${target} SYSTEM PRIVATE ${SPERR_INSTALL_PATH}/include)
    target_include_directories(${target} SYSTEM PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/../third_party/SZ3/include
        ${CMAKE_BINARY_DIR}/sz3_build/include
        ${CMAKE_CURRENT_SOURCE_DIR}/../third_party/SZ2/sz/include)
    target_compile_options(${target} PRIVATE -Wall -Ofast -march=native -Wno-dangling-else)
    target_link_libraries(${target} PRIVATE mrcz::mrcz zfp SZ3 SZ mrcz_qoz_bridge
        ${SPERR_INSTALL_PATH}/lib64/libSPERR.a ${MFEM_INSTALL_PATH}/lib/libmfem.a)
endfunction()

add_mrcz_executable(mrcz_compress compress.cpp)
add_mrcz_executable(mrcz_decompress decompress.cpp)

install(TARGETS mrcz_compress mrcz_decompress DESTINATION bin)
