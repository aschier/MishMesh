#include "MishMesh/smoothing.h"

#include <MishMesh/TriMesh.h>
#include <MishMesh/PolyMesh.h>

/**
 * Apply laplacian smoothing to a mesh by moving all vertices to the average of their neighbor's positions.
 * @param[inout] mesh The mesh to smooth.
 * @param[in] iterations The number of iterations for repeated smoothing of the mesh.
 * @param[in] component The components to be smoothed, i.e. tangential, normal, or both.
 * @note This is a simple interface to OpenMesh::Smoother. Use the OpenMesh functions directly,
 * when you need more control over the smoothing process.
 */
template<typename MeshT>
void MishMesh::smooth_mesh(MeshT &mesh, const int iterations, const typename MishMesh::SmootherT<MeshT>::Component component) {
	SmootherT<MeshT> smoother(mesh);
	smoother.initialize(component, OpenMesh::Smoother::JacobiLaplaceSmootherT<MeshT>::Continuity::C0);
	smoother.smooth(iterations);
}

template void MishMesh::smooth_mesh(MishMesh::TriMesh &mesh, const int iterations, const MishMesh::SmootherT<MishMesh::TriMesh>::Component component);
template void MishMesh::smooth_mesh(MishMesh::PolyMesh &mesh, const int iterations, const MishMesh::SmootherT<MishMesh::PolyMesh>::Component component);