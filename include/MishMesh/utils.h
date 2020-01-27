#pragma once
#include <MishMesh/TriMesh.h>
#include <MishMesh/PolyMesh.h>
#include <array>
#include <limits>

namespace MishMesh {
	TriMesh::HalfedgeHandle opposite_halfedge(const TriMesh &mesh, const TriMesh::FaceHandle &fh, const TriMesh::VertexHandle &vh);

	std::array<TriMesh::VertexHandle, 2> edge_vertices(const TriMesh &mesh, const TriMesh::HalfedgeHandle heh);
	std::array<TriMesh::VertexHandle, 2> edge_vertices(const TriMesh &mesh, const TriMesh::EdgeHandle eh);

	std::array<OpenMesh::Vec3d, 2> edge_points(const TriMesh &mesh, const TriMesh::HalfedgeHandle heh);
	std::array<OpenMesh::Vec3d, 2>  edge_points(const TriMesh &mesh, const TriMesh::EdgeHandle eh);

	std::array<TriMesh::VertexHandle, 3> face_vertices(const TriMesh &mesh, const TriMesh::FaceHandle fh);
	std::vector<PolyMesh::VertexHandle> face_vertices(const PolyMesh &mesh, const PolyMesh::FaceHandle fh);

	std::array<OpenMesh::Vec3d, 3> face_points(const TriMesh &mesh, const TriMesh::FaceHandle fh);
	std::vector<OpenMesh::Vec3d> face_points(const PolyMesh &mesh, const PolyMesh::FaceHandle fh);

	template<int DIM>
	double compute_area(const std::array<OpenMesh::VectorT<double, DIM>, 3> points);
	double compute_area(const TriMesh &mesh, const std::array<TriMesh::VertexHandle, 3> vertices);
	double compute_area(const TriMesh &mesh, const TriMesh::FaceHandle fh);

	TriMesh::VertexHandle obtuse_vertex(const TriMesh &mesh, const TriMesh::FaceHandle fh);


	template<typename MeshT>
	inline size_t euler_characteristic(const MeshT &mesh) {
		return mesh.n_vertices() - mesh.n_edges() + mesh.n_faces();
	}
}
