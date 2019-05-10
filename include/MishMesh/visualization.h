#pragma once
#include <MishMesh/TriMesh.h>
#include <vector>

namespace MishMesh{
	TriMesh edge_mesh(const TriMesh &mesh, std::vector<TriMesh::EdgeHandle> edge_handles, double thickness);
	TriMesh vertex_mesh(const TriMesh &mesh, std::vector<TriMesh::VertexHandle> vertex_handles, double size);
}