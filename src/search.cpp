#include "MishMesh/search.h"
#include <deque>

using namespace std;
using namespace MishMesh;

/**
 * Find faces connected to the given face using a breadth-first search.
 * @param input_mesh The mesh.
 * @param start_face The start face.
 * @returns A set of all faces reachable from the input face and the input face itself.
 */
set<TriMesh::FaceHandle> MishMesh::get_connected_faces(const TriMesh &input_mesh, const TriMesh::FaceHandle start_face){
	set<TriMesh::FaceHandle> faces{};
	list<TriMesh::FaceHandle> queue{start_face};
	while(!queue.empty()){
		auto face = queue.front();
		queue.pop_front();
		faces.insert(face);
		for(auto f_it = input_mesh.cff_begin(face); f_it != input_mesh.cff_end(face); f_it++){
			if(faces.find(*f_it) != faces.end()) continue;
			queue.push_back(*f_it);
		}
	}
	return faces;
}

/**
 * Find vertices connected to the given vertex using a breadth-first search.
 * @param input_mesh The mesh.
 * @param start_vertex The start vertex.
 * @returns A set of all faces reachable from the input vertex and the input vertex itself.
 */
set<TriMesh::VertexHandle> MishMesh::get_connected_vertices(const TriMesh &input_mesh, const TriMesh::VertexHandle start_vertex){
	set<TriMesh::VertexHandle> vertices{start_vertex};
	deque<TriMesh::VertexHandle> queue{start_vertex};
	while(!queue.empty()){
		auto vertex = queue.front();
		queue.pop_front();
		for(auto v_it = input_mesh.cvv_begin(vertex); v_it != input_mesh.cvv_end(vertex); v_it++){
			if(vertices.find(*v_it) != vertices.end()) continue;
			vertices.insert(*v_it);
			queue.push_back(*v_it);
		}
	}
	return vertices;
}

/**
 * Get the vertices of each connected component in the mesh.
 * @param input_mesh The mesh.
 * @returns A vector of sets of VertexHandles for the vertices of the connected components in the mesh.
 * @note The function does not guarantee a consistent order of the different connected groups with regard to
 * other functions dealing with connected components. Split the mesh and use get_connected_vertices on each submesh
 * when you need an ordered list of submeshes.
 */
std::vector<std::set<TriMesh::VertexHandle>> MishMesh::get_connected_components_vertices(const TriMesh &input_mesh) {
	std::vector<std::set<TriMesh::VertexHandle>> result;
	set<TriMesh::VertexHandle> vertices;
	for(auto v : input_mesh.vertices()) {
		vertices.insert(v);
	}
	while(!vertices.empty()) {
		const auto start_vertex = *vertices.begin();
		const auto component_vertices = get_connected_vertices(input_mesh, start_vertex);
		result.push_back(component_vertices);
		for(auto &vh : component_vertices) {
			vertices.erase(vh);
		}
	}
	return result;
}
