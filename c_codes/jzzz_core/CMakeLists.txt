add_library(libjzzz STATIC
    get_full_field_name.cpp
    read_block_info.cpp
    read_data_range.cpp
    read_header_file.cpp
)

# 设置 C++ 标准（与根项目一致）
set_target_properties(libjzzz PROPERTIES
    CXX_STANDARD 11
    CXX_STANDARD_REQUIRED ON
)

# 暴露头文件路径：让链接者（如主程序）能 #include "compressor.h"
target_include_directories(libjzzz PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# 可选：安装规则（make install 时，库去 lib/，头文件去 include/）
install(TARGETS libjzzz
    DESTINATION lib
)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} 
    DESTINATION include 
    FILES_MATCHING PATTERN "*.hpp"
)