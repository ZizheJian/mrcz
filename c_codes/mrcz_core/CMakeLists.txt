include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

add_library(mrcz INTERFACE)
add_library(mrcz::mrcz ALIAS mrcz)

target_include_directories(mrcz INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../third_party/get_max_cpus/c_codes>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_compile_features(mrcz INTERFACE cxx_std_20)

target_link_libraries(mrcz INTERFACE
    $<BUILD_INTERFACE:libzstd_static>
    $<INSTALL_INTERFACE:zstd::libzstd_static>
    tinyxml2::tinyxml2 mfem OpenMP::OpenMP_CXX)

install(TARGETS mrcz EXPORT MRCZTargets)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} 
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING PATTERN "*.hpp" PATTERN "*.tpp"
    PATTERN "parse_compress_args.hpp" EXCLUDE
    PATTERN "parse_decompress_args.hpp" EXCLUDE)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/../../third_party/get_max_cpus/c_codes/get_max_cpus.hpp DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

set(MRCZ_CMAKE_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/MRCZ")
install(EXPORT MRCZTargets FILE MRCZTargets.cmake NAMESPACE mrcz:: DESTINATION ${MRCZ_CMAKE_INSTALL_DIR})
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/../cmake/MRCZConfig.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/MRCZConfig.cmake INSTALL_DESTINATION ${MRCZ_CMAKE_INSTALL_DIR})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/MRCZConfig.cmake DESTINATION ${MRCZ_CMAKE_INSTALL_DIR})
