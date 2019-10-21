#include <OpenMesh/Core/IO/MeshIO.hh>
#include <MishMesh/TriMesh.h>
#include <MishMesh/sampling.h>
#include <OpenMesh/Core/Geometry/VectorT.hh>
#include <array>

int main(int argc, char **argv) {
	MishMesh::TriMesh mesh;
	std::array<OpenMesh::Vec3d, 3> triangle_points{
		OpenMesh::Vec3d{0.0, 0.0, 0.0},
		OpenMesh::Vec3d{1.0, 0.0, 0.0},
		OpenMesh::Vec3d{0.0, 1.0, 0.0}
	};
	for(int i = 0; i < 500; i++) {
		auto p = MishMesh::uniform_random_triangle_point(triangle_points);
		mesh.add_vertex(p);
	}
	OpenMesh::IO::write_mesh(mesh, "uniform_triangle.obj");

	mesh.clear();
	for(int i = 0; i < 500; i++) {
		auto p = MishMesh::uniform_random_annulus_points({0, 0}, 1.0, 0.5);
		mesh.add_vertex({p[0], p[1], 0.0});
	}
	OpenMesh::IO::write_mesh(mesh, "uniform_annulus.obj");
}
