cmake_minimum_required(VERSION 3.0.0)
project(MishMesh)
set(PROJECT_VERSION 0.1.3)

set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
include(CMakePackageConfigHelpers)
include(GNUInstallDirs)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

set(CMAKE_DEBUG_POSTFIX "-debug")
set(CMAKE_RELEASE_POSTFIX "-release")

option(MISHMESH_BUILD_EXAMPLE_PROGRAMS "Build MishMesh example programs" ON)
# When you include MishMesh as subdirectory you can prevent it from being
# installed together with your program by switching off its install targets.
option(MISHMESH_INSTALL_TARGETS "Create install targets for MishMesh." ON)

# Version 8.1 is required, but OpenMesh does not provide a package version file
find_package(OpenMesh REQUIRED)

find_package(Eigen3 3.3.7 QUIET)

find_package(OpenMP)
if (OPENMP_FOUND)
	set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_library(MishMesh
	include/MishMesh/OpenMesh/DoublePrecisionTraits.h
	include/MishMesh/TriMesh.h
	include/MishMesh/PolyMesh.h
	include/MishMesh/dijkstra.h
	src/dijkstra.cpp
	include/MishMesh/minimum_spanning_tree.h
	src/minimum_spanning_tree.cpp
	include/MishMesh/search.h
	src/search.cpp
	include/MishMesh/split.h
	src/split.cpp
	include/MishMesh/visualization.h
	src/visualization.cpp
	src/geodesics.cpp
	include/MishMesh/geodesics.h
	include/MishMesh/utils.h
	src/utils.cpp
	include/MishMesh/macros.h
	src/transformations.cpp
	include/MishMesh/transformations.h
	src/smoothing.cpp
	include/MishMesh/smoothing.h
	src/simplification.cpp
	include/MishMesh/simplification.h
	include/MishMesh/BBox.h
	include/MishMesh/sampling.h
	src/sampling.cpp
	include/MishMesh/poisson_disk_sampling.h
	src/poisson_disk_sampling.cpp
	include/MishMesh/VTK.h
	src/VTK.cpp
	include/MishMesh/OBJ.h
	src/OBJ.cpp
	include/MishMesh/dual_contouring.h
	src/meshing/dual_contouring.cpp
	src/meshing/vertexfit.hpp
	src/meshing/grid.hpp
)

target_compile_definitions(MishMesh PUBLIC -D_USE_MATH_DEFINES)
option(OM_STATIC_BUILD "Turn on when using a static OpenMesh libary" ON)
if(${OM_STATIC_BUILD})
	target_compile_definitions(MishMesh PUBLIC -DOM_STATIC_BUILD)
endif()

target_link_libraries(MishMesh PUBLIC
	OpenMeshCore
)
target_include_directories(MishMesh PUBLIC
	$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

if(${EIGEN3_FOUND})
	target_sources(MishMesh PRIVATE
		src/laplace.cpp
		src/heatGeodesics.cpp
		include/MishMesh/laplace.h
	)
	target_sources(MishMesh PRIVATE
		src/cone_singularities.cpp
		src/cone_singularities_impl.h
		include/MishMesh/cone_singularities.h
	)
	target_compile_definitions(MishMesh PUBLIC HAS_EIGEN)
	target_include_directories(MishMesh PUBLIC
		${EIGEN3_INCLUDE_DIR}
	)
else()
	message("Eigen3 was not found. Functions that need Eigen will not be available.")
endif()

if(MISHMESH_BUILD_EXAMPLE_PROGRAMS)
	# programs
	add_executable(split_connected_components
		src/programs/split_connected_components.cpp
		thirdparty/ProgramOptions.hxx
	)
	target_link_libraries(split_connected_components
		MishMesh)
	set_target_properties(split_connected_components PROPERTIES FOLDER programs)

	add_executable(geodesic_distances
		src/programs/geodesic_distances.cpp
		thirdparty/ProgramOptions.hxx
	)
	target_link_libraries(geodesic_distances
		MishMesh)
	set_target_properties(geodesic_distances PROPERTIES FOLDER programs)

	add_executable(mesh_smoother
		src/programs/mesh_smoother.cpp
		thirdparty/ProgramOptions.hxx
	)
	target_link_libraries(mesh_smoother
		MishMesh)

	set_target_properties(mesh_smoother PROPERTIES FOLDER programs)

	add_executable(grid
		src/programs/grid.cpp
	)
	target_link_libraries(grid
		MishMesh)
	set_target_properties(grid PROPERTIES FOLDER programs)

	add_executable(sampling_test
		src/programs/sampling_test.cpp
	)
	target_link_libraries(sampling_test
		MishMesh)
	set_target_properties(sampling_test PROPERTIES FOLDER programs)

	add_executable(meshing_test
		src/programs/meshing_test.cpp
	)
	target_link_libraries(meshing_test
		MishMesh)
	set_target_properties(meshing_test PROPERTIES FOLDER programs)

	add_executable(test_flatness
		src/programs/test_flatness.cpp
	)
	target_link_libraries(test_flatness
		MishMesh)
	set_target_properties(test_flatness PROPERTIES FOLDER programs)

	add_executable(icosphere
		src/programs/icosphere.cpp
	)
	target_link_libraries(icosphere
		MishMesh)
	set_target_properties(icosphere PROPERTIES FOLDER programs)

	if(${EIGEN3_FOUND})
		# Diffusion (laplace operator) example program
		add_executable(diffusion
			src/programs/diffusion.cpp
		)
		target_link_libraries(diffusion
			MishMesh)
		set_target_properties(diffusion PROPERTIES FOLDER programs)
		target_compile_definitions(diffusion PRIVATE HAS_EIGEN)

		# Compile definitions for geodesic_distances with heat geodesics
		target_compile_definitions(geodesic_distances PRIVATE HAS_EIGEN)

		# Cone Singularities example programs
		add_executable(cone_singularities
			src/programs/cone_singularities.cpp
		)
		target_link_libraries(cone_singularities
			MishMesh)
		set_target_properties(cone_singularities PROPERTIES FOLDER programs)
		add_executable(cone_cuts
			src/programs/cone_cuts.cpp
		)
		target_link_libraries(cone_cuts
			MishMesh)
		set_target_properties(cone_cuts PROPERTIES FOLDER programs)

		# install targets
		if(MISHMESH_INSTALL_TARGETS)
			install(TARGETS diffusion cone_singularities cone_cuts RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
		endif(MISHMESH_INSTALL_TARGETS)
	endif(${EIGEN3_FOUND})
endif(MISHMESH_BUILD_EXAMPLE_PROGRAMS)

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
	add_custom_target(doxygen ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} EXCLUDE_FROM_ALL)
	set_target_properties(doxygen PROPERTIES FOLDER "docs")
	add_custom_target(create_doc_directories COMMAND ${CMAKE_COMMAND} -E make_directory ${DOXYGEN_DIR} EXCLUDE_FROM_ALL)
	set_target_properties(create_doc_directories PROPERTIES FOLDER "docs")
	add_dependencies(doxygen create_doc_directories)
endif(DOXYGEN_FOUND)

if(MISHMESH_INSTALL_TARGETS)
	configure_package_config_file(MishMeshConfig.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/MishMeshConfig.cmake INSTALL_DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/MishMesh/cmake)
	write_basic_package_version_file( ${CMAKE_CURRENT_BINARY_DIR}/MishMeshConfigVersion.cmake COMPATIBILITY SameMajorVersion)
	install(TARGETS MishMesh EXPORT MishMesh-targets ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

	if(MISHMESH_BUILD_EXAMPLE_PROGRAMS)
		install(TARGETS split_connected_components mesh_smoother grid geodesic_distances RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
	endif(MISHMESH_BUILD_EXAMPLE_PROGRAMS)

	install(EXPORT MishMesh-targets FILE MishMeshTargets.cmake DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/MishMesh/cmake)
	install(FILES ${CMAKE_CURRENT_BINARY_DIR}/MishMeshConfig.cmake
		          ${CMAKE_CURRENT_BINARY_DIR}/MishMeshConfigVersion.cmake DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/MishMesh/cmake)
	install(DIRECTORY include/MishMesh DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} FILES_MATCHING PATTERN "*.h" PATTERN ".git" EXCLUDE)
endif(MISHMESH_INSTALL_TARGETS)
