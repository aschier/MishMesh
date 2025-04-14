#pragma once
#include "MishMesh/TriMesh.h"
#include "MishMesh/PolyMesh.h"

namespace MishMesh {
	void flip_with_uvs(MishMesh::TriMesh & mesh, const MishMesh::TriMesh::EdgeHandle &eh);
	void triangulate_with_uvs(MishMesh::PolyMesh &mesh);
}
