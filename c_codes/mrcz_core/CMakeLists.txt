add_library(mrcz INTERFACE)
add_library(mrcz::mrcz ALIAS mrcz)

# 暴露header-only MRCZ核心的头文件路径
target_include_directories(mrcz INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

# 优化选项
target_compile_options(mrcz INTERFACE -Wall -Ofast -march=native -fopenmp -Wno-dangling-else)

find_package(OpenMP REQUIRED)
target_link_libraries(mrcz INTERFACE libzstd_static tinyxml2::tinyxml2 OpenMP::OpenMP_CXX)

# INTERFACE target不生成.a或.so文件；安装时只登记target并复制公共头文件
install(TARGETS mrcz)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} 
    DESTINATION include 
    FILES_MATCHING PATTERN "*.hpp" PATTERN "*.tpp"
)
