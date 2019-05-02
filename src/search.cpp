#include "MishMesh/search.h"

using namespace std;
using namespace MishMesh;

/**
 * Find faces connected to the given face using a breadth-first search.
 * @param input_mesh The mesh.
 * @param start_face The start face.
 * @returns A set of all faces reachable from the input face and the input face itself.
 */
set<TriMesh::FaceHandle> MishMesh::find_connected_faces(const TriMesh &input_mesh, const TriMesh::FaceHandle start_face){
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
