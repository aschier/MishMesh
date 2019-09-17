#include "MishMesh/split.h"

#include <MishMesh/search.h>
#include <MishMesh/TriMesh.h>
#include <MishMesh/PolyMesh.h>
#include <MishMesh/macros.h>

using namespace std;
using namespace MishMesh;

/**
 * Build a mesh from a subset of faces of another mesh.
 * @param input_mesh The original mesh.
 * @param face_set A set of faces that should be included in the new mesh.
 * @param add_original_index_property When set to true, a property "orig_index" is added to the vertices and faces of the mesh.
 * @returns A MeshT mesh containing only the subset of the faces.
 * @note The method retains the connectivity of the original mesh,
 * but does not transfer any OpenMesh properties or other attributes.
 */
template<typename MeshT>
MeshT MishMesh::build_submesh(const MeshT &mesh, const set<typename MeshT::FaceHandle> &face_set, bool add_original_index_property) {
	MeshT submesh;
	OpenMesh::FPropHandleT<int> prop_orig_face_idx;
	OpenMesh::VPropHandleT<int> prop_orig_vertex_idx;
	if(add_original_index_property) {
		submesh.add_property(prop_orig_face_idx, "orig_index");
		submesh.add_property(prop_orig_vertex_idx, "orig_index");
	}
	map<typename MeshT::VertexHandle, typename MeshT::VertexHandle> vertex_map;
	for(auto &fh: face_set) {
		vector<typename MeshT::VertexHandle> vhs;
		FOR_CFV(v_it, fh) {
			if(vertex_map.find(*v_it) == vertex_map.end()) {
				const auto vh = submesh.add_vertex(mesh.point(*v_it));
				vertex_map[*v_it] = vh;
				if(add_original_index_property) {
					submesh.property(prop_orig_vertex_idx, vh) = v_it->idx();
				}
			}
			vhs.push_back(vertex_map[*v_it]);
		}
		const auto new_fh = submesh.add_face(vhs);
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
template<typename MeshT>
std::vector<MeshT> MishMesh::split_connected_components(const MeshT &input_mesh, bool add_original_index_property) {
	vector<MeshT> result_meshes;
	set<typename MeshT::FaceHandle> faces;
	for(auto f : input_mesh.faces()) {
		faces.insert(f);
	}
	while(!faces.empty()) {
		auto component_faces = get_connected_faces(input_mesh, *faces.begin());
		for(auto &f : component_faces){
			faces.erase(f);
		}
		MeshT component_mesh;
		component_mesh = build_submesh(input_mesh, component_faces, add_original_index_property);
		result_meshes.push_back(component_mesh);
	}
	return result_meshes;
}

template std::vector<MishMesh::TriMesh> MishMesh::split_connected_components(const MishMesh::TriMesh &input_mesh, bool add_original_index_property);
template std::vector<MishMesh::PolyMesh> MishMesh::split_connected_components(const MishMesh::PolyMesh &input_mesh, bool add_original_index_property);

template MishMesh::TriMesh MishMesh::build_submesh(const MishMesh::TriMesh &input_mesh, const set<typename MishMesh::TriMesh::FaceHandle> &face_set, bool add_original_index_property);
template MishMesh::PolyMesh MishMesh::build_submesh(const MishMesh::PolyMesh &input_mesh, const set<typename MishMesh::PolyMesh::FaceHandle> &face_set, bool add_original_index_property);
