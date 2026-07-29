cmake_minimum_required(VERSION 3.10)
project(mrcz)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build shared libraries" FORCE)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_subdirectory(mrcz_core)
add_subdirectory(../third_party/zstd ${CMAKE_BINARY_DIR}/zstd_build)
add_subdirectory(../third_party/zfp ${CMAKE_BINARY_DIR}/zfp_build)

set(SPERR_INSTALL_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../install/sperr_install")
set(MFEM_INSTALL_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../install/mfem_install")

function(add_mrcz_executable target source)
    add_executable(${target} ${source})
    target_include_directories(${target} PRIVATE ${SPERR_INSTALL_PATH}/include ${MFEM_INSTALL_PATH}/include)
    target_link_libraries(${target} PRIVATE libmrcz libzstd_static zfp
        ${SPERR_INSTALL_PATH}/lib64/libSPERR.a ${MFEM_INSTALL_PATH}/lib/libmfem.a)
endfunction()

add_mrcz_executable(mrcz_compress compress.cpp)
add_mrcz_executable(mrcz_decompress decompress.cpp)

install(TARGETS mrcz_compress mrcz_decompress DESTINATION bin)