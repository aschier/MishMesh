cmake_minimum_required(VERSION 3.30.0)
project(libacc)
include(CMakePackageConfigHelpers)
include(GNUInstallDirs)
set(CMAKE_CXX_STANDARD 17)

add_library(libacc INTERFACE)
target_include_directories(libacc INTERFACE
	$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

configure_package_config_file(libaccConfig.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/libaccConfig.cmake INSTALL_DESTINATION ${CMAKE_INSTALL_PREFIX})
install(TARGETS libacc EXPORT libacc-targets)

install(EXPORT libacc-targets FILE libaccTargets.cmake DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/libacc/cmake)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libaccConfig.cmake DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/libacc/cmake)
install(DIRECTORY "acc" DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} FILES_MATCHING PATTERN "*.h")
