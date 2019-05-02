#pragma once

#include <MishMesh/TriMesh.h>
#include <set>

namespace MishMesh {
	MishMesh::TriMesh build_submesh(const MishMesh::TriMesh &input_mesh, const std::set<MishMesh::TriMesh::FaceHandle>& face_set, bool add_original_index_property = false);
	std::vector<TriMesh> split_connected_components(const TriMesh &input_mesh, bool add_original_index_property = false);
}
