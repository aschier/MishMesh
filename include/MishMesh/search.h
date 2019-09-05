#pragma once

#include <set>
#include <list>
#include <vector>

namespace MishMesh {
	template<typename MeshT>
	std::set<typename MeshT::FaceHandle> get_connected_faces(const MeshT &input_mesh, const typename MeshT::FaceHandle start_face);
	template<typename MeshT>
	std::set<typename MeshT::VertexHandle> get_connected_vertices(const MeshT &input_mesh, const typename MeshT::VertexHandle start_vertex);
	template<typename MeshT>
	std::vector<std::set<typename MeshT::VertexHandle>> get_connected_components_vertices(const typename MeshT &input_mesh);
}
