#include "MishMesh/search.h"

#include <MishMesh/macros.h>
#include <MishMesh/TriMesh.h>
#include <MishMesh/PolyMesh.h>

#include <deque>

using namespace std;
using namespace MishMesh;

/**
 * Find faces connected to the given face using a breadth-first search.
 * @param mesh The mesh.
 * @param start_face The start face.
 * @returns A set of all faces reachable from the input face and the input face itself.
 */
template<typename MeshT>
set<typename MeshT::FaceHandle> MishMesh::get_connected_faces(const MeshT &mesh, const typename MeshT::FaceHandle start_face){
	set<typename MeshT::FaceHandle> faces{start_face};
	deque<typename MeshT::FaceHandle> queue{start_face};
	while(!queue.empty()){
		auto fh = queue.front();
		queue.pop_front();
		FOR_CFF(f_it, fh) {
			if(faces.find(*f_it) != faces.end()) continue;
			faces.insert(*f_it);
			queue.push_back(*f_it);
		}
	}
	return faces;
}

/**
 * Find vertices connected to the given vertex using a breadth-first search.
 * @param mesh The mesh.
 * @param start_vertex The start vertex.
 * @returns A set of all faces reachable from the input vertex and the input vertex itself.
 */
template<typename MeshT>
set<typename MeshT::VertexHandle> MishMesh::get_connected_vertices(const MeshT &mesh, const typename MeshT::VertexHandle start_vertex){
	set<typename MeshT::VertexHandle> vertices{start_vertex};
	deque<typename MeshT::VertexHandle> queue{start_vertex};
	while(!queue.empty()){
		auto vh = queue.front();
		queue.pop_front();
		FOR_CVV(v_it, vh) {
			if(vertices.find(*v_it) != vertices.end()) continue;
			vertices.insert(*v_it);
			queue.push_back(*v_it);
		}
	}
	return vertices;
}

/**
 * Get the vertices of each connected component in the mesh.
 * @param mesh The mesh.
 * @returns A vector of sets of VertexHandles for the vertices of the connected components in the mesh.
 * @note The function does not guarantee a consistent order of the different connected groups with regard to
 * other functions dealing with connected components. Split the mesh and use get_connected_vertices on each submesh
 * when you need an ordered list of submeshes.
 */
template<typename MeshT>
std::vector<std::set<typename MeshT::VertexHandle>> MishMesh::get_connected_components_vertices(const MeshT &mesh) {
	std::vector<std::set<typename MeshT::VertexHandle>> result;
	set<typename MeshT::VertexHandle> vertices;
	for(auto v : mesh.vertices()) {
		vertices.insert(v);
	}
	while(!vertices.empty()) {
		const auto start_vertex = *vertices.begin();
		const auto component_vertices = get_connected_vertices(mesh, start_vertex);
		result.push_back(component_vertices);
		for(auto &vh : component_vertices) {
			vertices.erase(vh);
		}
	}
	return result;
}

template set<typename TriMesh::FaceHandle> MishMesh::get_connected_faces(const TriMesh &mesh, const TriMesh::FaceHandle start_face);
template set<PolyMesh::FaceHandle> MishMesh::get_connected_faces(const PolyMesh &mesh, const PolyMesh::FaceHandle start_face);
template set<TriMesh::VertexHandle> MishMesh::get_connected_vertices(const TriMesh &mesh, const TriMesh::VertexHandle start_vertex);
template set<PolyMesh::VertexHandle> MishMesh::get_connected_vertices(const PolyMesh &mesh, const PolyMesh::VertexHandle start_vertex);
template std::vector<std::set<TriMesh::VertexHandle>> MishMesh::get_connected_components_vertices(const TriMesh &mesh);
template std::vector<std::set<PolyMesh::VertexHandle>> MishMesh::get_connected_components_vertices(const PolyMesh &mesh);
