#pragma once
#include "MishMesh/TriMesh.h"
#include "MishMesh/PolyMesh.h"
#include "MishMesh/utils.h"

namespace MishMesh {
	void flip_with_uvs(MishMesh::TriMesh &mesh, const MishMesh::TriMesh::EdgeHandle &eh);
	SmartHalfedgePair split_with_uvs(MishMesh::TriMesh &mesh, const OpenMesh::SmartHalfedgeHandle heh, const double t);
	void triangulate_with_uvs(MishMesh::PolyMesh &mesh);
}
