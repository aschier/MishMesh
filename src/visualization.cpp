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
				vhs[cube_indices[k][0]],
				vhs[cube_indices[k][1]],
				vhs[cube_indices[k][2]],
				vhs[cube_indices[k][3]]
			};
			auto fh = vertexMesh.add_face(face_vec);
		}
	}
	return vertexMesh;
}
