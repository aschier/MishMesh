MishMesh
========

A mishmash of useful mesh functions.

Build
-----
Set CMake variables:
- ``OPENMESH_LIBRARY_DIR``: The root folder of an OpenMesh installation.
- ``EIGEN3_DIR``: The CMake directory of a installed libeigen, e.g. ``<EIGEN_PREFIX>/share/eigen3/cmake``.

Optional Dependencies
---------------------
- When doxygen is found, a ``doxygen`` target is created to generate API documentation.
- When OpenMP is available, some algorithms are parallelized.
