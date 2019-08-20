#include "MishMesh/simplification.h"

#include <MishMesh/TriMesh.h>
#include <MishMesh/PolyMesh.h>

/**
 * Collapse short edges to improve the triangle quality.
 */
template<typename MeshT>
void MishMesh::collapse_short_edges(MeshT &mesh, const double epsilon) {
	mesh.request_edge_status();
	mesh.request_face_status();
	mesh.request_vertex_status();
	for(auto eh : mesh.edges()) {
		auto heh = mesh.halfedge_handle(eh, 0);
		if(mesh.calc_edge_length(heh) < epsilon){
			mesh.collapse(heh);
		}
	}
	mesh.garbage_collection();
	mesh.release_edge_status();
	mesh.release_face_status();
	mesh.release_vertex_status();
}

template void MishMesh::collapse_short_edges(TriMesh &mesh, const double epsilon);
template void MishMesh::collapse_short_edges(PolyMesh &mesh, const double epsilon);