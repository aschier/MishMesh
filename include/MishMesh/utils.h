#pragma once
#include <MishMesh/TriMesh.h>
#include <array>
#include <limits>

#include <MishMesh/BBox.h>

namespace MishMesh {
	TriMesh::HalfedgeHandle opposite_halfedge(const TriMesh &mesh, const TriMesh::FaceHandle &fh, TriMesh::VertexHandle &vh);
	std::array<TriMesh::VertexHandle, 3> face_vertices(const TriMesh &mesh, const TriMesh::FaceHandle fh);

	template<int DIM>
	double compute_area(const std::array<OpenMesh::VectorT<double, DIM>, 3> points);
	double compute_area(const TriMesh &mesh, const std::array<TriMesh::VertexHandle, 3> vertices);
	double compute_area(const TriMesh &mesh, const TriMesh::FaceHandle fh);
}
