#include "MishMesh/transformations.h"
#include "MishMesh/PolyMesh.h"
#include "MishMesh/utils.h"

using namespace std;

/**
 * Embed the triangle in the XY plane, with p_1=(0,0), p_2=(0, |p3D_2 - p3D_1|)
 *
 * @param p1 The first point.
 * @param p2 The second point.
 * @param p3 The third point.
 * @returns an array of 2D points.
 */
array<OpenMesh::Vec2d, 3> MishMesh::embed_triangle(OpenMesh::Vec3d p1, OpenMesh::Vec3d p2, OpenMesh::Vec3d p3) {
    OpenMesh::Vec3d x = (p2 - p1).normalized();
    OpenMesh::Vec3d n = (x % (p3 - p1)).normalized();
    OpenMesh::Vec3d y = n % x;
    return std::array<OpenMesh::Vec2d, 3> {
        OpenMesh::Vec2d {0.0, 0.0},
        OpenMesh::Vec2d {(p2 - p1).norm(), 0.0},
        OpenMesh::Vec2d {x.dot(p3 - p1), y.dot(p3 - p1)}
    };
}

/**
 * Embed the triangle in the XY plane, with p_1=(0,0), p_2=(0, |p3D_2 - p3D_1|)
 *
 * @param p1 The first point.
 * @param p2 The second point.
 * @param p3 The third point.
 * @returns an array of 2D points.
 */
array<OpenMesh::Vec2d, 3> MishMesh::embed_triangle(OpenMesh::Vec2d p1, OpenMesh::Vec2d p2, OpenMesh::Vec2d p3) {
    OpenMesh::Vec2d x = (p2 - p1).normalized();
    OpenMesh::Vec2d y = {-x[1], x[0]};
    return std::array<OpenMesh::Vec2d, 3> {
        OpenMesh::Vec2d {0.0, 0.0},
        OpenMesh::Vec2d {(p2 - p1).norm(), 0.0},
        OpenMesh::Vec2d {x.dot(p3 - p1), y.dot(p3 - p1)}
    };
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

std::array<OpenMesh::Vec2d, 3> MishMesh::embed_triangle(TriMesh &mesh, const TriMesh::HalfedgeHandle heh) {
	OpenMesh::Vec3d p = mesh.point(mesh.opposite_vh(heh));
	return embed_triangle(mesh, heh, p);
}


template std::array<OpenMesh::Vec2d, 3> MishMesh::embed_triangle(std::array<OpenMesh::VectorT<double, 2>, 3> points);
template std::array<OpenMesh::Vec2d, 3> MishMesh::embed_triangle(std::array<OpenMesh::VectorT<double, 3>, 3> points);
template std::array<OpenMesh::Vec2d, 3> MishMesh::embed_triangle(MishMesh::TriMesh &mesh, const typename MishMesh::TriMesh::HalfedgeHandle heh, const OpenMesh::Vec3d p);
template std::array<OpenMesh::Vec2d, 3> MishMesh::embed_triangle(MishMesh::PolyMesh &mesh, const typename MishMesh::PolyMesh::HalfedgeHandle heh, const OpenMesh::Vec3d p);
