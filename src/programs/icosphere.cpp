#include <MishMesh/TriMesh.h>
#include <MishMesh/PolyMesh.h>
#include <MishMesh/macros.h>
#include <MishMesh/utils.h>
#include <MishMesh/OBJ.h>

#include <vector>
#include <array>
#include <iostream>

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

/**
 * Global in-place 1 to 4 mesh refinement
 * @param mesh The mesh
 * @param iterations Number of refinement iterations
 */
void refine_1_to_4(MishMesh::TriMesh &mesh, int iterations = 1) {
	OpenMesh::EPropHandleT<OpenMesh::VertexHandle> prop_new_vh;
	mesh.add_property(prop_new_vh);
	mesh.request_face_status();
	mesh.request_edge_status();
	mesh.request_vertex_status();
	for(int iteration = 0; iteration < iterations; iteration++) {
		for(auto eh : mesh.edges()) {
			auto p = mesh.calc_edge_midpoint(eh);
			mesh.property(prop_new_vh, eh) = mesh.add_vertex(p / p.norm());
		}
		std::vector<std::array<OpenMesh::VertexHandle, 3>> new_faces;
		for(auto fh : mesh.faces()) {
			FOR_CFH(h_it, fh) {
				new_faces.push_back({mesh.property(prop_new_vh, h_it->edge()),
				                     h_it->to(),
				                     mesh.property(prop_new_vh, h_it->next().edge())});
			}
			new_faces.push_back({mesh.property(prop_new_vh, fh.halfedge().edge()),
			                     mesh.property(prop_new_vh, fh.halfedge().next().edge()),
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

MishMesh::PolyMesh dualsphere(MishMesh::TriMesh &icosphere) {
	MishMesh::PolyMesh mesh;
	OpenMesh::FPropHandleT<OpenMesh::VertexHandle> prop_dual_vh;
	icosphere.add_property(prop_dual_vh);

	for(auto fh : icosphere.faces()) {
		const std::array<MishMesh::TriMesh::VertexHandle, 3> vhs = MishMesh::face_vertices(icosphere, fh);
		auto p = (icosphere.point(vhs[0]) + icosphere.point(vhs[1]) + icosphere.point(vhs[2])) / 3.0;
		p /= p.norm(); // normalize the sphere to radius 1
		icosphere.property(prop_dual_vh, fh) = mesh.add_vertex(p);
	}
	for(auto vh : icosphere.vertices()) {
		std::vector<MishMesh::PolyMesh::VertexHandle> vhs;
		for(auto f_it = icosphere.cvf_ccwbegin(vh); f_it != icosphere.cvf_ccwend(vh); f_it++) {
			vhs.push_back(icosphere.property(prop_dual_vh, *f_it));
		}
		mesh.add_face(vhs.data(), vhs.size());
	}

	icosphere.remove_property(prop_dual_vh);
	return mesh;
}

int main(int argc, char **argv) {
	int iterations = 0;
	if(argc > 1) {
		iterations = atoi(argv[1]);
	}
	MishMesh::TriMesh mesh = create_isosahedron();
	refine_1_to_4(mesh, iterations);
	std::cerr << mesh.n_vertices() << " " << mesh.n_edges() << " " << mesh.n_faces() << std::endl;
	mesh.release_face_normals();
	MishMesh::writeOBJ(mesh, "icosphere.obj");
	auto dualmesh = dualsphere(mesh);
	dualmesh.release_face_normals();
	MishMesh::writeOBJ(dualmesh, "dualsphere.obj");
}
