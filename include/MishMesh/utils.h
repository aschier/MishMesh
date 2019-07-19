#pragma once
#include <MishMesh/TriMesh.h>
#include <array>

namespace MishMesh {
	TriMesh::HalfedgeHandle opposite_halfedge(const TriMesh &mesh, const TriMesh::FaceHandle &fh, TriMesh::VertexHandle &vh);
	std::array<TriMesh::VertexHandle, 3> face_vertices(const TriMesh &mesh, const TriMesh::FaceHandle fh);
}