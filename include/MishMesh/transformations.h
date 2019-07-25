#pragma once

#include <MishMesh/TriMesh.h>
#include <array>

namespace MishMesh {
	template<int DIM>
	std::array<OpenMesh::Vec2d, 3> embed_triangle(OpenMesh::VectorT<double, DIM> p1, OpenMesh::VectorT<double, DIM> p2, OpenMesh::VectorT<double, DIM> p3);

	template std::array<OpenMesh::Vec2d, 3> embed_triangle(OpenMesh::VectorT<double, 2> p1, OpenMesh::VectorT<double, 2> p2, OpenMesh::VectorT<double, 2> p3);
	template std::array<OpenMesh::Vec2d, 3> embed_triangle(OpenMesh::VectorT<double, 3> p1, OpenMesh::VectorT<double, 3> p2, OpenMesh::VectorT<double, 3> p3);
}