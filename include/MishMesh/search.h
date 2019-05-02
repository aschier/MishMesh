#pragma once

#include <MishMesh/TriMesh.h>
#include <set>
#include <list>

namespace MishMesh {
	std::set<MishMesh::TriMesh::FaceHandle> find_connected_faces(const MishMesh::TriMesh &input_mesh, const MishMesh::TriMesh::FaceHandle start_face);
}
