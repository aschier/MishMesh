#pragma once

#include <MishMesh/TriMesh.h>
#include <array>

#ifdef HAS_EIGEN
#include <Eigen/Eigen>
#endif

namespace MishMesh {
	std::array<OpenMesh::Vec2d, 3> embed_triangle(OpenMesh::Vec3d p1, OpenMesh::Vec3d p2, OpenMesh::Vec3d p3);
	std::array<OpenMesh::Vec2d, 3> embed_triangle(OpenMesh::Vec2d p1, OpenMesh::Vec2d p2, OpenMesh::Vec2d p3);
	template<int DIM>
	std::array<OpenMesh::Vec2d, 3> embed_triangle(std::array<OpenMesh::VectorT<double, DIM>, 3> points);
	template<typename MeshT>
	std::array<OpenMesh::Vec2d, 3> embed_triangle(MeshT &mesh, const typename MeshT::HalfedgeHandle heh, const OpenMesh::Vec3d p);
	std::array<OpenMesh::Vec2d, 3> embed_triangle(TriMesh &mesh, const TriMesh::HalfedgeHandle heh);

#ifdef HAS_EIGEN
	std::array<Eigen::Vector2d, 3> embed_triangle(const Eigen::Vector3d &p1, const Eigen::Vector3d &p2, const Eigen::Vector3d &p3);
	std::array<Eigen::Vector2d, 3> embed_triangle(const std::array<Eigen::Vector3d, 3> &points);
#endif

}