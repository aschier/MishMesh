#include "MishMesh/transformations.h"
#include "MishMesh/PolyMesh.h"

using namespace std;

/**
 * Embed the triangle in the XY plane, with p_1=(0,0), p_2=(0, |p3D_2 - p3D_1|) and p_3 projected into the positive half plane.
 *
 * @param p1 The first point.
 * @param p2 The second point.
 * @param p3 The third point.
 * @returns an array of 2D points.
 */
template<int DIM>
array<OpenMesh::Vec2d, 3> MishMesh::embed_triangle(OpenMesh::VectorT<double, DIM> p1, OpenMesh::VectorT<double, DIM> p2, OpenMesh::VectorT<double, DIM> p3) {
	double v2x = (p2 - p1).norm();
	assert(v2x != 0);
	double v3x = (p2 - p1) | (p3 - p1) / (p2 - p1).norm();
	double C = pow((p3 - p1).norm(), 2);
	double A = pow(v3x, 2);
	double v3y = C > A ? sqrt(C-A) : sqrt(A-C);
	OpenMesh::Vec2d v1 = {0, 0};
	OpenMesh::Vec2d v2 = {v2x, 0};
	OpenMesh::Vec2d v3 = {v3x, v3y};
	return {v1, v2, v3};
}

template<int DIM>
std::array<OpenMesh::Vec2d, 3> MishMesh::embed_triangle(std::array<OpenMesh::VectorT<double, DIM>, 3> points) {
	return embed_triangle(points[0], points[1], points[2]);
}

template<typename MeshT>
std::array<OpenMesh::Vec2d, 3> MishMesh::embed_triangle(MeshT &mesh, const typename MeshT::HalfedgeHandle heh, const OpenMesh::Vec3d p) {
	auto p1 = mesh.point(mesh.from_vertex_handle(heh));
	auto p2 = mesh.point(mesh.to_vertex_handle(heh));
	return MishMesh::embed_triangle(p1, p2, p);
}

template std::array<OpenMesh::Vec2d, 3> MishMesh::embed_triangle(OpenMesh::VectorT<double, 2> p1, OpenMesh::VectorT<double, 2> p2, OpenMesh::VectorT<double, 2> p3);
template std::array<OpenMesh::Vec2d, 3> MishMesh::embed_triangle(OpenMesh::VectorT<double, 3> p1, OpenMesh::VectorT<double, 3> p2, OpenMesh::VectorT<double, 3> p3);
template std::array<OpenMesh::Vec2d, 3> MishMesh::embed_triangle(std::array<OpenMesh::VectorT<double, 2>, 3> points);
template std::array<OpenMesh::Vec2d, 3> MishMesh::embed_triangle(std::array<OpenMesh::VectorT<double, 3>, 3> points);
template std::array<OpenMesh::Vec2d, 3> MishMesh::embed_triangle(MishMesh::TriMesh &mesh, const typename MishMesh::TriMesh::HalfedgeHandle heh, const OpenMesh::Vec3d p);
template std::array<OpenMesh::Vec2d, 3> MishMesh::embed_triangle(MishMesh::PolyMesh &mesh, const typename MishMesh::PolyMesh::HalfedgeHandle heh, const OpenMesh::Vec3d p);
