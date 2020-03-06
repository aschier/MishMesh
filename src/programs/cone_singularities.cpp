#include <iostream>
#include <set>

#include <OpenMesh/Core/IO/MeshIO.hh>

#include <MishMesh/cone_singularities.h>
#include <MishMesh/visualization.h>
#include <MishMesh/BBox.h>

using namespace std;

int main(int argc, char **argv) {
	MishMesh::TriMesh mesh;
	std::string filename(argv[1]);
	OpenMesh::IO::read_mesh(mesh, filename);

	MishMesh::BBox<> bbox = MishMesh::bounding_box(mesh);
	const double diameter = bbox.diameter();

	auto singularity_vhs = MishMesh::compute_cone_singularities(mesh, 1.0, 10);

	auto singularity_mesh = MishMesh::vertex_mesh(mesh, singularity_vhs, diameter / 100.0);
	OpenMesh::IO::write_mesh(singularity_mesh, filename + "_singularities.obj");
}
