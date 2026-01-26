cmake_minimum_required(VERSION 3.10)
project(jzzz)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build shared libraries" FORCE)

# LTO 修复：设置 ar/ranlib（在 flags 前）
# find_program(CMAKE_AR NAMES gcc-ar gar)
# if (NOT CMAKE_AR)
#     message(FATAL_ERROR "gcc-ar not found. Install gcc-ar.")
# endif()
# find_program(CMAKE_RANLIB NAMES gcc-ranlib ranlib)
# if (NOT CMAKE_RANLIB)
#     message(FATAL_ERROR "gcc-ranlib not found.")
# endif()

# 优化标志
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Ofast -march=native -flto")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Ofast -march=native -fopenmp -mavx2 -mfma")

add_subdirectory(jzzz_core)
add_subdirectory(third_party/zstd)
add_subdirectory(third_party/zfp)

add_executable(jzzz_compress compress.cpp)
target_include_directories(jzzz_compress PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/third_party/sperr/install/include)
# target_link_directories(jzzz_compress PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/third_party/sperr/install/lib64)
target_link_libraries(jzzz_compress PRIVATE libjzzz libzstd_static zfp ${CMAKE_CURRENT_SOURCE_DIR}/third_party/sperr/install/lib64/libSPERR.a)
    
add_executable(jzzz_decompress decompress.cpp)
target_include_directories(jzzz_decompress PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/third_party/sperr/install/include)
# target_link_directories(jzzz_decompress PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/third_party/sperr/install/lib64)
target_link_libraries(jzzz_decompress PRIVATE libjzzz libzstd_static zfp ${CMAKE_CURRENT_SOURCE_DIR}/third_party/sperr/install/lib64/libSPERR.a)

install(TARGETS jzzz_compress jzzz_decompress DESTINATION bin)