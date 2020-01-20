#pragma once

#include <MishMesh/TriMesh.h>
#include <array>

namespace MishMesh {
	template<int DIM>
	std::array<OpenMesh::Vec2d, 3> embed_triangle(OpenMesh::VectorT<double, DIM> p1, OpenMesh::VectorT<double, DIM> p2, OpenMesh::VectorT<double, DIM> p3);
	template<int DIM>
	std::array<OpenMesh::Vec2d, 3> embed_triangle(std::array<OpenMesh::VectorT<double, DIM>, 3> points);
	template<typename MeshT>
	std::array<OpenMesh::Vec2d, 3> embed_triangle(MeshT &mesh, const typename MeshT::HalfedgeHandle heh, const OpenMesh::Vec3d p);
}