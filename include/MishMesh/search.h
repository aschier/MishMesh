#pragma once

#include <MishMesh/TriMesh.h>
#include <set>
#include <list>

namespace MishMesh {
	std::set<MishMesh::TriMesh::FaceHandle> get_connected_faces(const MishMesh::TriMesh &input_mesh, const MishMesh::TriMesh::FaceHandle start_face);
	std::set<TriMesh::VertexHandle> get_connected_vertices(const MishMesh::TriMesh &input_mesh, const MishMesh::TriMesh::VertexHandle start_vertex);
	std::vector<std::set<TriMesh::VertexHandle>> get_connected_components_vertices(const TriMesh &input_mesh);
}
