include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

add_library(mrcz INTERFACE)
add_library(mrcz::mrcz ALIAS mrcz)

target_include_directories(mrcz INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_compile_features(mrcz INTERFACE cxx_std_20)

target_link_libraries(mrcz INTERFACE
    $<BUILD_INTERFACE:libzstd_static>
    $<INSTALL_INTERFACE:zstd::libzstd_static>
    get_max_cpus::get_max_cpus tinyxml2::tinyxml2 mfem OpenMP::OpenMP_CXX)

# get_max_cpus的头文件由它自己的CMakeLists装到include/get_max_cpus，这里把它的目标一起导出，外部项目链接mrcz::mrcz时才能找到这个头文件
install(TARGETS mrcz libgetmaxcpus EXPORT MRCZTargets)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING PATTERN "*.hpp" PATTERN "*.tpp"
    PATTERN "parse_compress_args.hpp" EXCLUDE
    PATTERN "parse_decompress_args.hpp" EXCLUDE)

set(MRCZ_CMAKE_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/MRCZ")
install(EXPORT MRCZTargets FILE MRCZTargets.cmake NAMESPACE mrcz:: DESTINATION ${MRCZ_CMAKE_INSTALL_DIR})
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/../cmake/MRCZConfig.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/MRCZConfig.cmake INSTALL_DESTINATION ${MRCZ_CMAKE_INSTALL_DIR})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/MRCZConfig.cmake DESTINATION ${MRCZ_CMAKE_INSTALL_DIR})
