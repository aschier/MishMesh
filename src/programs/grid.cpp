#include <cstdlib>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <MishMesh/TriMesh.h>
#include <MishMesh/visualization.h>
#include <MishMesh/BBox.h>

int main(int argc, char **argv) {
	int resolution[3] = {10, 10, 10};
	if(argc > 1) {
		int res = std::atoi(argv[1]);
		resolution[0] = res;
		resolution[1] = res;
		resolution[2] = res;
	}
	MishMesh::BBox<> bbox{{0,0,0}, {1,1,1}};
	auto mesh = MishMesh::grid_mesh(resolution, bbox, 0.2 / resolution[0]);
	OpenMesh::IO::write_mesh(mesh, "grid.obj");
}
