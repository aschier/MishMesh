#pragma once

#include <MishMesh/TriMesh.h>

namespace MishMesh {
	typedef OpenMesh::VPropHandleT<double> GeodesicDistanceProperty;
	void compute_geodesics(TriMesh &mesh, const TriMesh::VertexHandle start_vh, const GeodesicDistanceProperty geodesicDistanceProperty);
}