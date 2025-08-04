#pragma once
#include <MishMesh/TriMesh.h>
#include <MishMesh/BBox.h>
#include <vector>

namespace MishMesh {
	template<typename MeshT>
	void add_box(MeshT &mesh, const BBox<OpenMesh::Vec3d, 3> box);
	TriMesh edge_mesh(const TriMesh &mesh, const std::vector<TriMesh::EdgeHandle> &edge_handles, const double thickness, const OpenMesh::EPropHandleT<double> *prop_edge_value = nullptr);
	TriMesh edge_mesh(const TriMesh &mesh, const OpenMesh::EPropHandleT<bool> prop_add_edge, const double thickness, const OpenMesh::EPropHandleT<double> *prop_edge_value = nullptr);
	TriMesh edge_mesh(TriMesh &mesh, const std::vector<TriMesh::EdgeHandle> &edge_handles, const double thickness, const OpenMesh::EPropHandleT<double> *prop_edge_value, const double offset = 0.0);
	TriMesh vertex_mesh(const TriMesh &mesh, std::vector<TriMesh::VertexHandle> vertex_handles, double size);
	template<typename MeshT>
	void colorize_mesh(MeshT &mesh, const OpenMesh::VPropHandleT<double> &vertexProperty);
	template<typename MeshT>
	void colorize_mesh(MeshT &mesh, const OpenMesh::FPropHandleT<double> &faceProperty);
	void cosine_colorize_mesh(MishMesh::TriMesh &mesh, const OpenMesh::VPropHandleT<double> &vertexProperty, const double periods);
	TriMesh grid_mesh(const int resolution[3], const BBox<OpenMesh::Vec3d, 3> bbox, const double point_size = 0);
	TriMesh isosurface_grid_mesh(const int resolution[3], const BBox<OpenMesh::Vec3d, 3> bbox, std::vector<double> &point_values, double point_size);
}
