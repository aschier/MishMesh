#include <MishMesh/visualization.h>
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

constexpr short cube_indices[6][4]{
	{0, 1, 3, 2},
	{1, 0, 4, 5},
	{0, 2, 6, 4},
	{2, 3, 7, 6},
	{3, 1, 5, 7},
	{5, 4, 6, 7},
};

/**
 * Create a mesh, that visualizes a set of vertices from an input mesh as boxes.
 * @param mesh The input mesh.
 * @param vertex_handles The vertices to visualize.
 * @param size The side length of the boxes.
 */
TriMesh MishMesh::vertex_mesh(const TriMesh &mesh, std::vector<TriMesh::VertexHandle> vertex_handles, double size) {
	TriMesh vertexMesh;
	for(auto vh : vertex_handles) {
		array<TriMesh::VertexHandle, 8> vhs;
		for(short j = 0; j < 8; j++) {
			auto p = mesh.point(vh);
			vhs[j] = vertexMesh.add_vertex({
				p[0] + ((j & 1) == 1 ? size : -size) / 2.0,
				p[1] + ((j >> 1 & 1) == 1 ? size : -size) / 2.0,
				p[2] + ((j >> 2 & 1) == 1 ? size : -size) / 2.0
				});
		}
		for(short k = 0; k < 6; k++) {
			vector<TriMesh::VertexHandle> face_vec{
				vhs[cube_indices[k][3]],
				vhs[cube_indices[k][2]],
				vhs[cube_indices[k][1]],
				vhs[cube_indices[k][0]]
			};
			auto fh = vertexMesh.add_face(face_vec);
		}
	}
	return vertexMesh;
}

/**
 * Colorize the mesh from black to red using a given vertex property. The colors will be scaled from minimum to maximum value.
 * @param[inout] mesh The mesh.
 * @param[in] vertexProperty The vertex property.
 * @note You need to request_vertex_color before using this method.
 */
void MishMesh::colorize_mesh(MishMesh::TriMesh &mesh, const OpenMesh::VPropHandleT<double> &vertexProperty) {
	assert(mesh.has_vertex_colors());
	double max_value = -numeric_limits<double>::infinity();
	double min_value = numeric_limits<double>::infinity();
	for(auto vh : mesh.vertices()) {
		double value = mesh.property(vertexProperty, vh);
		max_value = std::max(max_value, value);
		min_value = std::min(min_value, value);
	}
	for(auto vh : mesh.vertices()) {
		double value = mesh.property(vertexProperty, vh);
		mesh.set_color(vh, {255 * value / (max_value - min_value), 0, 0});
	}
}

/**
 * Generate a regular grid of a bounding box defined by its top-left-far and right-bottom-near vertices.
 * @param resolution The number of grid points in each dimension.
 * @param bbox_ltf The left-top-far vertex of the bounding box.
 * @param bbox_rbn The right-bottom-near vertex of the bounding box.
 * @param point_size The size of a cube at a grid point. Use 0 to use 1/10 of the width of a grid cell.
 * @returns A vertex mesh with cubes at the grid points.
 */
TriMesh MishMesh::grid_mesh(const int resolution[3], const double bbox_ltf[3], const double bbox_rbn[3], double point_size) {
	MishMesh::TriMesh point_mesh;
	vector<MishMesh::TriMesh::VertexHandle> vertices;
	vertices.reserve(resolution[0] * resolution[1] * resolution[2]);
	for(int i = 0; i < resolution[0]; i++) {
		for(int j = 0; j < resolution[1]; j++) {
			for(int k = 0; k < resolution[2]; k++){
				OpenMesh::Vec3d point{
					bbox_ltf[0] + (bbox_rbn[0] - bbox_ltf[0]) / resolution[0] * i,
					bbox_ltf[1] + (bbox_rbn[1] - bbox_ltf[1]) / resolution[0] * j,
					bbox_ltf[2] + (bbox_rbn[2] - bbox_ltf[2]) / resolution[0] * k,
				};
				vertices.push_back(point_mesh.add_vertex(point));
			}
		}
	}
	if(point_size == 0){
		point_size = abs(bbox_rbn[0] - bbox_ltf[0]) / resolution[0] / 10.0;
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
 * @param bbox_ltf The left-top-far vertex of the bounding box.
 * @param bbox_rbn The right-bottom-near vertex of the bounding box.
 * @param point_values The function values on the grid vertices, using the index format x + y*resolution[0] + z*resolution[0]*resolution[1].
 * @param point_size The size of a cube at a grid point. Use 0 to use 1/10 of the width of a grid cell.
 * @returns A vertex mesh with cubes at the grid points.
 */
TriMesh MishMesh::isosurface_grid_mesh(const int resolution[3], const double bbox_ltf[3], const double bbox_rbn[3], std::vector<double> &point_values, double point_size) {
	MishMesh::TriMesh point_mesh;
	vector<MishMesh::TriMesh::VertexHandle> vertices;
	for(int i = 1; i < resolution[0] - 1; i++) {
		for(int j = 1; j < resolution[1] - 1; j++) {
			for(int k = 1; k < resolution[2] - 1; k++){
				if(
					point_values[grid_cell_index(resolution, i + 1, j, k)] * point_values[grid_cell_index(resolution, i, j, k)] <= 0 ||
					point_values[grid_cell_index(resolution, i - 1, j, k)] * point_values[grid_cell_index(resolution, i, j, k)] <= 0 ||
					point_values[grid_cell_index(resolution, i, j + 1, k)] * point_values[grid_cell_index(resolution, i, j, k)] <= 0 ||
					point_values[grid_cell_index(resolution, i, j - 1, k)] * point_values[grid_cell_index(resolution, i, j, k)] <= 0 ||
					point_values[grid_cell_index(resolution, i, j, k + 1)] * point_values[grid_cell_index(resolution, i, j, k)] <= 0 ||
					point_values[grid_cell_index(resolution, i, j, k - 1)] * point_values[grid_cell_index(resolution, i, j, k)] <= 0
					) {
					OpenMesh::Vec3d p = OpenMesh::Vec3d{
						bbox_ltf[0] + i * (bbox_rbn[0] - bbox_ltf[0]) / resolution[0],
						bbox_ltf[1] + j * (bbox_rbn[1] - bbox_ltf[1]) / resolution[1],
						bbox_ltf[2] + k * (bbox_rbn[2] - bbox_ltf[2]) / resolution[2]};
					vertices.push_back(point_mesh.add_vertex(p));
				}
			}
		}
	}
	if(point_size == 0){
		point_size = abs(bbox_rbn[0] - bbox_ltf[0]) / resolution[0] / 10.0;
	}
	return MishMesh::vertex_mesh(point_mesh, vertices, point_size);
}