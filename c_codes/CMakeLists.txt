cmake_minimum_required(VERSION 3.10)
project(jzzz)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

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
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3 -march=native -fopenmp -mavx2 -mfma")

add_subdirectory(jzzz_core)  # 你的子目录
add_subdirectory(third_party/zstd)

add_executable(jzzz_compress compress.cpp)
target_link_libraries(jzzz_compress libjzzz libzstd_static)
add_executable(jzzz_decompress decompress.cpp)
target_link_libraries(jzzz_decompress libjzzz libzstd_static)

install(TARGETS jzzz_compress jzzz_decompress DESTINATION bin)