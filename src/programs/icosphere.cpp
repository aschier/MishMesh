#include <MishMesh/TriMesh.h>
#include <MishMesh/PolyMesh.h>
#include <MishMesh/macros.h>
#include <MishMesh/utils.h>
#include <MishMesh/OBJ.h>
#include <MishMesh/refine.h>

#include <vector>
#include <array>
#include <iostream>

#include "dualsphere.hpp"

/**
 * Create a simple icosahedron mesh
 * @returns A icosahedron mesh
 */
MishMesh::TriMesh create_isosahedron() {
	const double PHI = (1.0 + sqrt(5)) / 2.0;
	MishMesh::TriMesh mesh;
	std::array<OpenMesh::VertexHandle, 12> vhs{
	    mesh.add_vertex({-1, PHI, 0}),
	    mesh.add_vertex({1, PHI, 0}),
	    mesh.add_vertex({-1, -PHI, 0}),
	    mesh.add_vertex({1, -PHI, 0}),

	    mesh.add_vertex({0, -1, PHI}),
	    mesh.add_vertex({0, 1, PHI}),
	    mesh.add_vertex({0, -1, -PHI}),
	    mesh.add_vertex({0, 1, -PHI}),

	    mesh.add_vertex({PHI, 0, -1}),
	    mesh.add_vertex({PHI, 0, 1}),
	    mesh.add_vertex({-PHI, 0, -1}),
	    mesh.add_vertex({-PHI, 0, 1})};

	mesh.add_face(vhs[0], vhs[11], vhs[5]);
	mesh.add_face(vhs[0], vhs[5], vhs[1]);
	mesh.add_face(vhs[0], vhs[1], vhs[7]);
	mesh.add_face(vhs[0], vhs[7], vhs[10]);
	mesh.add_face(vhs[0], vhs[10], vhs[11]);

	mesh.add_face(vhs[1], vhs[5], vhs[9]);
	mesh.add_face(vhs[5], vhs[11], vhs[4]);
	mesh.add_face(vhs[11], vhs[10], vhs[2]);
	mesh.add_face(vhs[10], vhs[7], vhs[6]);
	mesh.add_face(vhs[7], vhs[1], vhs[8]);

	mesh.add_face(vhs[3], vhs[9], vhs[4]);
	mesh.add_face(vhs[3], vhs[4], vhs[2]);
	mesh.add_face(vhs[3], vhs[2], vhs[6]);
	mesh.add_face(vhs[3], vhs[6], vhs[8]);
	mesh.add_face(vhs[3], vhs[8], vhs[9]);

	mesh.add_face(vhs[4], vhs[9], vhs[5]);
	mesh.add_face(vhs[2], vhs[4], vhs[11]);
	mesh.add_face(vhs[6], vhs[2], vhs[10]);
	mesh.add_face(vhs[8], vhs[6], vhs[7]);
	mesh.add_face(vhs[9], vhs[8], vhs[1]);

	// Normalize to a sphere of radius 1
	for(auto vh : mesh.vertices()) {
		mesh.point(vh) /= mesh.point(vh).norm();
	}

	return mesh;
}

int main(int argc, char **argv) {
	int iterations = 0;
	if(argc > 1) {
		iterations = atoi(argv[1]);
	}
	MishMesh::TriMesh mesh = create_isosahedron();
	MishMesh::global_refine_1_to_4(mesh, iterations);
	std::cerr << mesh.n_vertices() << " " << mesh.n_edges() << " " << mesh.n_faces() << std::endl;
	mesh.release_face_normals();
	MishMesh::writeOBJ(mesh, "icosphere.obj");
	auto dualmesh = dualsphere(mesh);
	dualmesh.release_face_normals();
	MishMesh::writeOBJ(dualmesh, "dualsphere.obj");
}
