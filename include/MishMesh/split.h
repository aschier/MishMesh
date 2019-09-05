#pragma once

#include <set>
#include <vector>

namespace MishMesh {
	template<typename MeshT>
	MeshT build_submesh(const MeshT &input_mesh, const std::set<typename MeshT::FaceHandle> &face_set, bool add_original_index_property = false);
	template<typename MeshT>
	std::vector<MeshT> split_connected_components(const MeshT &input_mesh, bool add_original_index_property = false);
}