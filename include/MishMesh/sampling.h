#pragma once

#include <array>
#include <OpenMesh/Core/Geometry/VectorT.hh>
#include <MishMesh/TriMesh.h>

namespace MishMesh {
	OpenMesh::Vec3d uniform_random_barycentric_coordinate(const std::array<OpenMesh::Vec3d, 3>);
	OpenMesh::Vec3d uniform_random_barycentric_coordinate(const TriMesh &mesh, TriMesh::FaceHandle fh);
	OpenMesh::Vec3d uniform_random_triangle_point(const std::array<OpenMesh::Vec3d, 3>);
	OpenMesh::Vec3d uniform_random_triangle_point(const TriMesh &mesh, TriMesh::FaceHandle fh);
	OpenMesh::Vec2d uniform_random_circle_point(const OpenMesh::Vec2d center, const double radius);
	OpenMesh::Vec2d uniform_random_annulus_point(const OpenMesh::Vec2d center, const double outer_radius, const double inner_radius);
}
