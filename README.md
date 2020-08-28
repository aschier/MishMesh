MishMesh
========

A mishmash of useful mesh functions.

Dependencies
------------

- OpenMESH 8.1

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

----------------------------------------------------------------------

Functions
---------

##### Data Types
* DoublePrecisionTraits allow for meshes using double instead of float
* ``TriMesh`` and ``PolyMesh`` types allow for consistent mesh types across different projects linking against each other.

##### Searching
* An implementation of Dijkstra's algorithm (shortest edge path search), that allows to use custom distance functions
* Breadth first search for elements
  * Get faces reachable from a source face
  * Get vertices reachable from a source vertex
  * Get the sets of vertices of the connected components in a mesh.

##### Minimum Spanning Trees
* Calculate a Minimum Spanning Tree of the Dijkstra paths between a set of vertices
* Calculate one Dijkstra Minimum Spanning Tree per connected connected component
* Calculate all shortest distances for a single source vertex

##### Splitting
* Build a submesh from a set of faces of a mesh
* Split mesh into its connected components

##### Sampling
* Uniform sampling inside triangles
* Uniform sampling inside circles
* Uniform sampling inside annuli (rings)
* Poisson disk sampling inside triangles and circles using the algorithm of Bridson (``Bridson, R. (2007). Fast Poisson disk sampling in arbitrary dimensions.
	 * In ACM SIGGRAPH 2007 Sketches on - SIGGRAPH ’07, (San Diego, California: ACM Press), pp. 22-es.``) [[URL](https://www.cs.ubc.ca/~rbridson/docs/bridson-siggraph07-poissondisk.pdf)]
* Fibonnaci sphere sampling

##### Visualization
* Create an "edge mesh", that highlights a set of edges of a mesh
* Create a "vertex mesh", that highlights a set of vertices of a mesh
* Colorize a mesh using a given vertex ``double`` property
* Construct a grid
* Construct a grid around an isosurface of a given signed distance function

##### Simplification & Smoothing
* A simple wrapper around the OpenMesh smoothing function
* Collapse short edges to improve mesh quality

##### Utilities
* Get an array of the vertex handles of a face
* Compute the triangle area of a face or a set of points or vertex handles
* Get the halfedge opposite to a vertex in a given face
* Compute the euler characteristic of the mesh
* Embed a 3D triangle in the 2D plane

##### Macros
* Add macros for often used iterator loops

**Example:**
```
MishMesh::TriMesh::FaceHandle fh = mesh.face_handle(0);
FOR_CFV(v_it, fh) {
    std::cout << "Vertex Index: " << v_it->idx() << std::endl;
}
```

##### Bounding Box
* A templated bounding box implementation with functions for testing if points are inside the box and for clipping points at the bounding box

##### Mesh Laplace
* Compute a mesh laplacian with support for normalization and vertex area weighting

##### Geodesics
* Compute geodesic distances using the algorithm of Novotni and Klein (``Novotni, M., & Klein, R. (2002). Computing geodesic distances on triangular meshes. In In Proc. of WSCG’2002``) [[URL](http://cg.cs.uni-bonn.de/de/publikationen/paper-details/novotni-2002-computing/)].
* Compute geodesic distances using the heat method (``Crane, K., Weischedel, C., and Wardetzky, M. (2013). Geodesics in Heat. ACM Trans. Graph. 32, 1–11``) [[URL](https://dl.acm.org/citation.cfm?id=2516977)].

##### Cone Singularities
* Compute cone singularities on the mesh using the algorithm of Ben-Chen et. al (``Ben-Chen, M., Gotsman, C., & Bunin, G. (2008). Conformal Flattening by Curvature Prescription and Metric Scaling. Computer Graphics Forum, 27(2), 449–458``) [[URL](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1467-8659.2008.01142.x)].
