#include "MishMesh/split.h"

#include <MishMesh/search.h>

using namespace std;
using namespace MishMesh;

/**
 * Build a mesh from a subset of faces of another mesh.
 * @param input_mesh The original mesh.
 * @param face_set A set of faces that should be included in the new mesh.
 * @param add_original_index_property When set to true, a property "orig_index" is added to the vertices and faces of the mesh.
 * @returns A TriMesh containing only the subset of the faces.
 * @note The method retains the connectivity of the original mesh,
 * but does not transfer any OpenMesh properties or other attributes.
 */
TriMesh MishMesh::build_submesh(const TriMesh &input_mesh, const set<TriMesh::FaceHandle> &face_set, bool add_original_index_property) {
	TriMesh submesh;
	OpenMesh::FPropHandleT<int> prop_orig_face_idx;
	OpenMesh::VPropHandleT<int> prop_orig_vertex_idx;
	if(add_original_index_property) {
		submesh.add_property(prop_orig_face_idx, "orig_index");
		submesh.add_property(prop_orig_vertex_idx, "orig_index");
	}
	map<TriMesh::VertexHandle, TriMesh::VertexHandle> vertex_map;
	for(auto &fh: face_set) {
		array<TriMesh::VertexHandle, 3> vhs;
		short i = 0;
		for(auto v_it = input_mesh.cfv_ccwbegin(fh); v_it != input_mesh.cfv_ccwend(fh); v_it++) {
			assert(i < 3);
			if(vertex_map.find(*v_it) == vertex_map.end()) {
				const auto vh = submesh.add_vertex(input_mesh.point(*v_it));
				vertex_map[*v_it] = vh;
				if(add_original_index_property) {
					submesh.property(prop_orig_vertex_idx, vh) = v_it->idx();
				}
			}
			vhs[i++] = vertex_map[*v_it];
		}
		const auto new_fh = submesh.add_face(vhs.data(), 3);
		if(add_original_index_property) {
			submesh.property(prop_orig_face_idx, new_fh) = fh.idx();
		}
	}
	return submesh;
}

/**
 * Split the mesh into its connected components.
 * @param input_mesh the mesh.
 * @param add_original_index_property When set to true, a property "orig_index" is added to the vertices and faces of the mesh.
 * @returns A vector with one mesh per connected component.
 */
std::vector<TriMesh> MishMesh::split_connected_components(const TriMesh &input_mesh, bool add_original_index_property) {
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
		component_mesh = build_submesh(input_mesh, component_faces, add_original_index_property);
		result_meshes.push_back(component_mesh);
	}
	return result_meshes;
}