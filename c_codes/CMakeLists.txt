cmake_minimum_required(VERSION 3.10)
project(mrcz)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build shared libraries" FORCE)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# 优化标志
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-unused-variable")

add_subdirectory(mrcz_core)
add_subdirectory(third_party/zstd)
add_subdirectory(third_party/zfp)

set(SPERR_INSTALL_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../install/sperr_install")
set(MFEM_INSTALL_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../install/mfem_install")

add_executable(mrcz_compress compress.cpp)
target_include_directories(mrcz_compress PRIVATE ${SPERR_INSTALL_PATH}/include)
target_link_libraries(mrcz_compress PRIVATE libmrcz libzstd_static zfp ${SPERR_INSTALL_PATH}/lib64/libSPERR.a)
target_include_directories(mrcz_compress PRIVATE ${MFEM_INSTALL_PATH}/include)
target_link_libraries(mrcz_compress PRIVATE ${MFEM_INSTALL_PATH}/lib/libmfem.a)
target_compile_options(mrcz_compress PRIVATE -Wall -Ofast -march=native -fopenmp -mavx2 -mfma -Wno-dangling-else)
target_link_options(mrcz_compress PRIVATE -fopenmp)
    
add_executable(mrcz_decompress decompress.cpp)
target_include_directories(mrcz_decompress PRIVATE ${SPERR_INSTALL_PATH}/include)
target_link_libraries(mrcz_decompress PRIVATE libmrcz libzstd_static zfp ${SPERR_INSTALL_PATH}/lib64/libSPERR.a)
target_include_directories(mrcz_decompress PRIVATE ${MFEM_INSTALL_PATH}/include)
target_link_libraries(mrcz_decompress PRIVATE ${MFEM_INSTALL_PATH}/lib/libmfem.a)
target_compile_options(mrcz_decompress PRIVATE -Wall -Ofast -march=native -fopenmp -mavx2 -mfma -Wno-dangling-else)
target_link_options(mrcz_decompress PRIVATE -fopenmp)

install(TARGETS mrcz_compress mrcz_decompress DESTINATION bin)