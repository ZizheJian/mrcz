cmake_minimum_required(VERSION 3.10)
project(mrcz)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build shared libraries" FORCE)

# 优化标志
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Ofast -march=native -fopenmp -mavx2 -mfma -Wno-dangling-else")

add_subdirectory(mrcz_core)
add_subdirectory(third_party/zstd)
add_subdirectory(third_party/zfp)

add_executable(mrcz_compress compress.cpp)
target_include_directories(mrcz_compress PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/third_party/sperr/install/include)
# target_link_directories(mrcz_compress PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/third_party/sperr/install/lib64)
target_link_libraries(mrcz_compress PRIVATE libmrcz libzstd_static zfp ${CMAKE_CURRENT_SOURCE_DIR}/third_party/sperr/install/lib64/libSPERR.a)
    
add_executable(mrcz_decompress decompress.cpp)
target_include_directories(mrcz_decompress PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/third_party/sperr/install/include)
# target_link_directories(mrcz_decompress PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/third_party/sperr/install/lib64)
target_link_libraries(mrcz_decompress PRIVATE libmrcz libzstd_static zfp ${CMAKE_CURRENT_SOURCE_DIR}/third_party/sperr/install/lib64/libSPERR.a)

install(TARGETS mrcz_compress mrcz_decompress DESTINATION bin)