MishMesh
========

A mishmash of useful mesh functions.

Dependencies
------------

- OpenMESH

Optional Dependencies
---------------------

- doxygen: When doxygen is found, a ``doxygen`` target is created to generate API documentation.
- OpenMP: When OpenMP is available, some algorithms are parallelized.
- Eigen3: Needed for computing laplace matrices.

Build
-----

Set CMake variables:
- ``OPENMESH_LIBRARY_DIR``: The root folder of an OpenMesh installation.
- ``Eigen3_DIR``: The CMake directory of a installed libeigen, e.g. ``<EIGEN_PREFIX>/share/eigen3/cmake``.
