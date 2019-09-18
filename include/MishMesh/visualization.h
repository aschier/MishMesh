#pragma once
#include <MishMesh/TriMesh.h>
#include <MishMesh/BBox.h>
#include <vector>

namespace MishMesh{
	template<typename MeshT>
	void add_box(MeshT &mesh, const BBox<OpenMesh::Vec3d, 3> box);
	TriMesh edge_mesh(const TriMesh &mesh, std::vector<TriMesh::EdgeHandle> edge_handles, double thickness);
	TriMesh vertex_mesh(const TriMesh &mesh, std::vector<TriMesh::VertexHandle> vertex_handles, double size);
	void colorize_mesh(MishMesh::TriMesh &mesh, const OpenMesh::VPropHandleT<double> &vertexProperty);
	TriMesh grid_mesh(const int resolution[3], const BBox<OpenMesh::Vec3d, 3> bbox, const double point_size = 0);
	TriMesh isosurface_grid_mesh(const int resolution[3], const BBox<OpenMesh::Vec3d, 3> bbox, std::vector<double> &point_values, double point_size);
}
