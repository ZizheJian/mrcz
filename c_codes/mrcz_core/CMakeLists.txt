add_library(libmrcz INTERFACE)

# 暴露头文件路径：让链接者（如主程序）能 #include "compressor.h"
target_include_directories(libmrcz INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

# 优化选项
target_compile_options(libmrcz INTERFACE -Wall -Ofast -march=native -fopenmp -Wno-dangling-else)

find_package(OpenMP REQUIRED)
target_link_libraries(libmrcz INTERFACE libzstd_static tinyxml2::tinyxml2 OpenMP::OpenMP_CXX)

# 可选：安装规则（make install 时，库去 lib/，头文件去 include/）
install(TARGETS libmrcz DESTINATION lib)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} 
    DESTINATION include 
    FILES_MATCHING PATTERN "*.hpp" PATTERN "*.tpp"
)
