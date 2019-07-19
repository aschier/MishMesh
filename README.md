MishMesh
========

A mishmash of useful mesh functions.

Build
-----
Set CMake variables:
- ``OPENMESH_LIBRARY_DIR``: The root folder of an OpenMesh installation.

Optional Dependencies
---------------------
- When doxygen is found, a ``doxygen`` target is created to generate API documentation.
- When OpenMP is available, some algorithms are parallelized.
