cmake_minimum_required(VERSION 2.8.10)
project(MishMesh)
include(CMakePackageConfigHelpers)
include(GNUInstallDirs)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

set(CMAKE_DEBUG_POSTFIX "d")

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
find_package(OpenMesh 8.0   REQUIRED)
find_package(Eigen3   3.3.7 REQUIRED)

add_library(MishMesh
	include/MishMesh/OpenMesh/DoublePrecisionTraits.h
	include/MishMesh/TriMesh.h
	include/MishMesh/dijkstra.h
	include/MishMesh/minimum_spanning_tree.h
	src/dijkstra.cpp
	src/minimum_spanning_tree.cpp
)
target_compile_definitions(MishMesh PUBLIC -D_USE_MATH_DEFINES)
target_link_libraries(MishMesh PUBLIC
	${OPENMESH_LIBRARIES}
)
target_include_directories(MishMesh PUBLIC
	$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
	${OPENMESH_INCLUDE_DIR}
	${EIGEN3_INCLUDE_DIR}
)

# doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
	set(DOXYGEN_DIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DOCDIR})
	if(DOXYGEN_DOT_EXECUTABLE)
		set(HAVE_DOT "YES")
	else()
		set(HAVE_DOT "NO")
	endif()
	configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
	add_custom_target(doxygen ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
	set_target_properties(doxygen PROPERTIES FOLDER "docs")
	add_custom_target(create_doc_directories COMMAND ${CMAKE_COMMAND} -E make_directory ${DOXYGEN_DIR} EXCLUDE_FROM_ALL)
	set_target_properties(create_doc_directories PROPERTIES FOLDER "docs")
	add_dependencies(doxygen create_doc_directories)
endif(DOXYGEN_FOUND)

configure_package_config_file(MishMeshConfig.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/MishMeshConfig.cmake INSTALL_DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/MishMesh/cmake)
install(TARGETS MishMesh EXPORT MishMesh-targets ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

install(EXPORT MishMesh-targets FILE MishMeshTargets.cmake DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/MishMesh/cmake)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/MishMeshConfig.cmake DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/MishMesh/cmake)
install(DIRECTORY include/MishMesh DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} FILES_MATCHING PATTERN "*.h" PATTERN ".git" EXCLUDE)
