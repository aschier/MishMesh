#include "MishMesh/split.h"

#include <MishMesh/search.h>

using namespace std;
using namespace MishMesh;

/**
 * Build a mesh from a subset of faces of another mesh.
 * @param input_mesh The original mesh.
 * @param face_set A set of faces that should be included in the new mesh.
 * @returns A TriMesh containing only the subset of the faces.
 * @note The method retains the connectivity of the original mesh,
 * but does not transfer any OpenMesh properties or other attributes.
 */
TriMesh MishMesh::build_submesh(const TriMesh &input_mesh, const set<TriMesh::FaceHandle> &face_set) {
	TriMesh submesh;
	map<TriMesh::VertexHandle, TriMesh::VertexHandle> vertex_map;
	for(auto &f: face_set) {
		array<TriMesh::VertexHandle, 3> vhs;
		short i = 0;
		for(auto v_it = input_mesh.cfv_ccwbegin(f); v_it != input_mesh.cfv_ccwend(f); v_it++) {
			assert(i < 3);
			if(vertex_map.find(*v_it) == vertex_map.end()) {
				vertex_map[*v_it] = submesh.add_vertex(input_mesh.point(*v_it));
			}
			vhs[i++] = vertex_map[*v_it];
		}
		submesh.add_face(vhs.data(), 3);
	}
	return submesh;
}

/**
 * Split the mesh into its connected components.
 * @param input_mesh the mesh.
 * @returns A vector with one mesh per connected component.
 */
std::vector<TriMesh> MishMesh::split_connected_components(const TriMesh &input_mesh) {
	vector<TriMesh> result_meshes;
	set<TriMesh::FaceHandle> faces;
	for(auto f : input_mesh.faces()) {
		faces.insert(f);
	}
	while(!faces.empty()) {
		auto component_faces = find_connected_faces(input_mesh, *faces.begin());
		for(auto &f : component_faces){
			faces.erase(f);
		}
		TriMesh component_mesh;
		component_mesh = build_submesh(input_mesh, component_faces);
		result_meshes.push_back(component_mesh);
	}
	return result_meshes;
}