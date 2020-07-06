#include "MishMesh/visualization.h"
#include <MishMesh/PolyMesh.h>
#include <MishMesh/BBox.h>

#include <array>

using namespace MishMesh;
using namespace std;

/**
 * Create a mesh, that visualizes selected edges from an input mesh as rectangles.
 * @param mesh The input mesh.
 * @param edge_handles The edges to visualize.
 * @param thickness The relative thickness of the edges.
 * @returns The edge mesh.
 */
TriMesh MishMesh::edge_mesh(const TriMesh &mesh, std::vector<TriMesh::EdgeHandle> edge_handles, double thickness) {
	TriMesh edgeMesh;
	for(auto eh : edge_handles) {
		std::array<TriMesh::VertexHandle, 4> vhs;
		auto heh = mesh.halfedge_handle(eh, 0);
		auto vh1 = mesh.from_vertex_handle(heh);
		auto vh2 = mesh.to_vertex_handle(heh);
		OpenMesh::Vec3d edge = mesh.point(vh2) - mesh.point(vh1);
		OpenMesh::Vec3d d = OpenMesh::cross(edge, OpenMesh::Vec3d{0,0,1});
		if(d == OpenMesh::Vec3d{0, 0, 0}) {
			d = OpenMesh::cross(edge, OpenMesh::Vec3d{0,1,0});
		}
		d *= thickness / d.norm();
		vhs[0] = edgeMesh.add_vertex(mesh.point(vh1));
		vhs[1] = edgeMesh.add_vertex(mesh.point(vh2));
		vhs[2] = edgeMesh.add_vertex(mesh.point(vh1) + d);
		vhs[3] = edgeMesh.add_vertex(mesh.point(vh2) + d);

		auto fh1 = edgeMesh.add_face(vhs[0], vhs[1], vhs[2]);
		auto fh2 = edgeMesh.add_face(vhs[2], vhs[1], vhs[3]);
	}
	return edgeMesh;
}

constexpr short box_indices[6][4]{
	{0, 1, 3, 2},
	{1, 0, 4, 5},
	{0, 2, 6, 4},
	{2, 3, 7, 6},
	{3, 1, 5, 7},
	{5, 4, 6, 7},
};

/**
 * Add a box given by a MishMesh::Bbox object to a mesh.
 * @param mesh The mesh.
 * @param box The box.
 * @tparam MeshT the type of the mesh, e.g. MishMesh::TriMesh or MishMesh::PolyMesh.
 *         The faces are added as quads, but types like MishMesh::TriMesh triangulate them
 *         automatically.
 */
template<typename MeshT>
void MishMesh::add_box(MeshT &mesh, const BBox<OpenMesh::Vec3d, 3> box) {
		array<typename MeshT::VertexHandle, 8> vhs;
		for(short j = 0; j < 8; j++) {
			vhs[j] = mesh.add_vertex({
				((j & 1) == 0 ? box.ltf[0] : box.rbn[0]),
				((j >> 1 & 1) == 0 ? box.ltf[1] : box.rbn[1]),
				((j >> 2 & 1) == 0 ? box.ltf[2] : box.rbn[2]),
				});
		}
		for(short k = 0; k < 6; k++) {
			vector<typename MeshT::VertexHandle> face_vec{
				vhs[box_indices[k][3]],
				vhs[box_indices[k][2]],
				vhs[box_indices[k][1]],
				vhs[box_indices[k][0]]
			};
			mesh.add_face(face_vec);
		}
}

/**
 * Create a mesh, that visualizes a set of vertices from an input mesh as boxes.
 * @param mesh The input mesh.
 * @param vertex_handles The vertices to visualize.
 * @param size The side length of the boxes.
 */
TriMesh MishMesh::vertex_mesh(const TriMesh &mesh, std::vector<TriMesh::VertexHandle> vertex_handles, double size) {
	TriMesh vertexMesh;
	for(auto vh : vertex_handles) {
		auto p = mesh.point(vh);
		BBox<OpenMesh::Vec3d, 3> cube {
			p - OpenMesh::Vec3d{size, size, size},
			p + OpenMesh::Vec3d{size, size, size},
		};
		add_box(vertexMesh, cube);
	}
	return vertexMesh;
}

/**
 * Colorize the mesh from black to red using a given vertex property. The colors are be scaled from minimum to maximum value.
 *
 * @param[inout] mesh The mesh.
 * @param[in] vertexProperty The vertex property.
 * @note You need to use request_vertex_color before using this method.
 */
template<typename MeshT>
void MishMesh::colorize_mesh(MeshT &mesh, const OpenMesh::VPropHandleT<double> &vertexProperty) {
	assert(mesh.has_vertex_colors());
	double max_value = -numeric_limits<double>::infinity();
	double min_value = numeric_limits<double>::infinity();
	for(auto vh : mesh.vertices()) {
		double value = mesh.property(vertexProperty, vh);
		if(!isfinite(value)) continue;
		max_value = std::max(max_value, value);
		min_value = std::min(min_value, value);
	}
	for(auto vh : mesh.vertices()) {
		double value = mesh.property(vertexProperty, vh);
		if(isfinite(value)) {
			mesh.set_color(vh, {255 * (value - min_value) / (max_value - min_value), 0, 0});
		} else {
			mesh.set_color(vh, {0, 0, 255});
		}
	}
}

/**
 * Colorize the mesh using a cosine function to create a repeating gradient red to black using a given vertex property.
 * The values are be scaled from minimum to maximum value, so the number of cosine periods between minimum and maximum
 * can be chosen using the periods parameter.
 *
 * @param[inout] mesh The mesh.
 * @param[in] vertexProperty The vertex property.
 * @param[in] periods The number of periods for the cosine function.
 * @note You need to use request_vertex_color before using this method.
 */
void MishMesh::cosine_colorize_mesh(MishMesh::TriMesh &mesh, const OpenMesh::VPropHandleT<double> &vertexProperty, const double periods) {
	assert(mesh.has_vertex_colors());
	double max_value = -numeric_limits<double>::infinity();
	double min_value = numeric_limits<double>::infinity();
	for(auto vh : mesh.vertices()) {
		double value = mesh.property(vertexProperty, vh);
		if(!isfinite(value)) continue;
		max_value = std::max(max_value, value);
		min_value = std::min(min_value, value);
	}
	for(auto vh : mesh.vertices()) {
		double value = ((mesh.property(vertexProperty, vh) - min_value) / (max_value - min_value)) * 2 * M_PI * periods;
		if(isfinite(value)) {
			mesh.set_color(vh, {static_cast<unsigned char>(255 * (1.0 + cos(value)) / 2.0), 0, 0});
		} else {
			mesh.set_color(vh, {0, 0, 255});
		}
	}
}

/**
 * Generate a regular grid of a bounding box defined by its top-left-far and right-bottom-near vertices.
 * @param resolution The number of grid points in each dimension.
 * @param bbox The bounding box.
 * @param point_size The size of a cube at a grid point. Use 0 to use 1/10 of the width of a grid cell.
 * @returns A vertex mesh with cubes at the grid points.
 */
TriMesh MishMesh::grid_mesh(const int resolution[3], const BBox<OpenMesh::Vec3d, 3> bbox, double point_size) {
	MishMesh::TriMesh point_mesh;
	vector<MishMesh::TriMesh::VertexHandle> vertices;
	vertices.reserve(resolution[0] * resolution[1] * resolution[2]);
	for(int i = 0; i < resolution[0]; i++) {
		for(int j = 0; j < resolution[1]; j++) {
			for(int k = 0; k < resolution[2]; k++){
				OpenMesh::Vec3d point{
					bbox.ltf[0] + (bbox.rbn[0] - bbox.ltf[0]) / resolution[0] * i,
					bbox.ltf[1] + (bbox.rbn[1] - bbox.ltf[1]) / resolution[0] * j,
					bbox.ltf[2] + (bbox.rbn[2] - bbox.ltf[2]) / resolution[0] * k,
				};
				vertices.push_back(point_mesh.add_vertex(point));
			}
		}
	}
	if(point_size == 0){
		point_size = abs(bbox.rbn[0] - bbox.ltf[0]) / resolution[0] / 10.0;
	}
	return MishMesh::vertex_mesh(point_mesh, vertices, point_size);
}

inline int grid_cell_index(const int resolution[3], const int i, const int j, const int k) {
	return i + resolution[0] * j + resolution[0] * resolution[1] * k;
}

/**
 * Visualize the grid near a isosurface of a signed distance function, by generating a regular grid and only
 * showing the vertices that have at least one neighbor for which the function value the opposite sign.
 * @param resolution The number of grid points in each dimension.
 * @param bbox The bounding box.
 * @param point_values The function values on the grid vertices, using the index format x + y*resolution[0] + z*resolution[0]*resolution[1].
 * @param point_size The size of a cube at a grid point. Use 0 to use 1/10 of the width of a grid cell.
 * @returns A vertex mesh with cubes at the grid points.
 */
TriMesh MishMesh::isosurface_grid_mesh(const int resolution[3], const BBox<OpenMesh::Vec3d, 3> bbox, std::vector<double> &point_values, double point_size) {
	MishMesh::TriMesh point_mesh;
	vector<MishMesh::TriMesh::VertexHandle> vertices;
	for(int i = 1; i < resolution[0] - 1; i++) {
		for(int j = 1; j < resolution[1] - 1; j++) {
			for(int k = 1; k < resolution[2] - 1; k++){
				double v000 = point_values[grid_cell_index(resolution, i, j, k)];
				if(!isfinite(v000)) {
					continue;
				}
				double v001 = point_values[grid_cell_index(resolution, i, j, k+1)];
				double v010 = point_values[grid_cell_index(resolution, i, j+1, k)];
				double v100 = point_values[grid_cell_index(resolution, i+1, j, k)];
				double v00_1 = point_values[grid_cell_index(resolution, i, j, k-1)];
				double v0_10 = point_values[grid_cell_index(resolution, i, j-1, k)];
				double v_100 = point_values[grid_cell_index(resolution, i-1, j, k)];
				if(
					(isfinite(v001) && v001 * v000 <= 0) ||
					(isfinite(v010) && v010 * v000 <= 0) ||
					(isfinite(v100) && v100 * v000 <= 0) ||
					(isfinite(v00_1) && v00_1 * v000 <= 0) ||
					(isfinite(v0_10) && v0_10 * v000 <= 0) ||
					(isfinite(v_100) && v_100 * v000 <= 0)
					) {
					OpenMesh::Vec3d p = OpenMesh::Vec3d{
						bbox.ltf[0] + i * (bbox.rbn[0] - bbox.ltf[0]) / resolution[0],
						bbox.ltf[1] + j * (bbox.rbn[1] - bbox.ltf[1]) / resolution[1],
						bbox.ltf[2] + k * (bbox.rbn[2] - bbox.ltf[2]) / resolution[2]};
					vertices.push_back(point_mesh.add_vertex(p));
				}
			}
		}
	}
	if(point_size == 0){
		point_size = abs(bbox.rbn[0] - bbox.ltf[0]) / resolution[0] / 10.0;
	}
	return MishMesh::vertex_mesh(point_mesh, vertices, point_size);
}

template void MishMesh::add_box(MishMesh::TriMesh &mesh, const BBox<OpenMesh::Vec3d, 3> box);
template void MishMesh::add_box(MishMesh::PolyMesh &mesh, const BBox<OpenMesh::Vec3d, 3> box);

template void MishMesh::colorize_mesh(MishMesh::TriMesh &mesh, const OpenMesh::VPropHandleT<double> &vertexProperty);
template void MishMesh::colorize_mesh(MishMesh::PolyMesh &mesh, const OpenMesh::VPropHandleT<double> &vertexProperty);