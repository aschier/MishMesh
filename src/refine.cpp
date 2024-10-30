#include "MishMesh/refine.h"
#include <MishMesh/macros.h>

/**
 * Global in-place 1 to 4 mesh refinement
 *
 * Refines the mesh faces by splitting the edge midpoints and creating four new faces inside each face.
 * This function operates on all faces simultaneously, preserving the angles in all faces.
 *
 * @param mesh The mesh to refine
 * @param iterations Number of refinement iterations
 */

void MishMesh::global_refine_1_to_4(MishMesh::TriMesh &mesh, int iterations) {
	OpenMesh::EPropHandleT<OpenMesh::VertexHandle> prop_new_vh;
	mesh.add_property(prop_new_vh);
	mesh.request_face_status();
	mesh.request_edge_status();
	mesh.request_vertex_status();

	for(int iteration = 0; iteration < iterations; ++iteration) {
		for(auto eh : mesh.edges()) {
			auto p = mesh.calc_edge_midpoint(eh);
			mesh.property(prop_new_vh, eh) = mesh.add_vertex(p / p.norm());
		}

		std::vector<std::array<OpenMesh::VertexHandle, 3>> new_faces;
		for(auto fh : mesh.faces()) {
			FOR_CFH(h_it, fh) {
				new_faces.push_back(
				    {mesh.property(prop_new_vh, h_it->edge()), h_it->to(), mesh.property(prop_new_vh, h_it->next().edge())});
			}
			new_faces.push_back({mesh.property(prop_new_vh, fh.halfedge().edge()), mesh.property(prop_new_vh, fh.halfedge().next().edge()),
			                     mesh.property(prop_new_vh, fh.halfedge().next().next().edge())});
		}

		for(auto fh : mesh.faces()) {
			mesh.delete_face(fh, false);
		}

		for(auto new_face : new_faces) {
			mesh.add_face(new_face.data(), 3);
		}
	}

	mesh.garbage_collection();
	mesh.remove_property(prop_new_vh);
	mesh.release_face_status();
	mesh.release_edge_status();
	mesh.release_vertex_status();
}
