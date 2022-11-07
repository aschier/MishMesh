#include <MishMesh/TriMesh.h>
#include <MishMesh/PolyMesh.h>
#include <MishMesh/macros.h>
#include <MishMesh/utils.h>
#include <MishMesh/OBJ.h>

#include <vector>
#include <array>
#include <iostream>

#include "dualsphere.hpp"

/**
 * Create a simple UV sphere
 * @note if u is smaller than 3 it is set to 3 and if v is smaller than 2 it is set to 2.
*/
MishMesh::PolyMesh uvsphere(uint u, uint v) {
	u = std::max<uint>(u, 3);
	v = std::max<uint>(v, 2);
	v = v-1;

    MishMesh::PolyMesh mesh;
    auto vh_top = mesh.add_vertex({0, 1, 0});
    auto vh_bottom = mesh.add_vertex({0, -1, 0});
    std::vector<MishMesh::PolyMesh::VertexHandle> vhs(u * v);
    for (uint j = 0; j < v; j++) {
        double theta = M_PI * static_cast<double>(j + 1) / (v + 1);
        for (uint i = 0; i < u; i++) {
            uint idx = j * u + i;
            double phi = M_PI * static_cast<double>(i) / u;
            vhs[idx] = mesh.add_vertex(
                {std::sin(theta) * std::cos(2 * phi), std::cos(theta), std::sin(theta) * std::sin(2 * phi)});
        }
    }

    // Caps
    for (uint i = 0; i < u; i++) {
        mesh.add_face(vhs[(i + 1) % u], vhs[(i + 0) % u], vh_top);
        mesh.add_face(vhs[u * (v - 1) + (i + 0) % u], vhs[u * (v - 1) + (i + 1) % u], vh_bottom);
    }
    // Sides
    for (int j = 0; j < v - 1; j++) {
        for (int i = 0; i < u; i++) {
            mesh.add_face(vhs[(j + 0) % v * u + (i + 0) % u], vhs[(j + 0) % v * u + (i + 1) % u],
                          vhs[(j + 1) % v * u + (i + 1) % u], vhs[(j + 1) % v * u + (i + 0) % u]);
        }
    }

    return mesh;
}

int main(int argc, char **argv) {
    uint u = 20;
    uint v = 20;
    if (argc > 2) {
        u = atoi(argv[1]);
        v = atoi(argv[2]);
    }
    MishMesh::PolyMesh mesh = uvsphere(u, v);
    std::cerr << mesh.n_vertices() << " " << mesh.n_edges() << " " << mesh.n_faces() << std::endl;
    mesh.release_face_normals();
    MishMesh::writeOBJ(mesh, "uvsphere.obj");
    auto dualmesh = dualsphere(mesh);
    dualmesh.release_face_normals();
    MishMesh::writeOBJ(dualmesh, "dualsphere.obj");
}
