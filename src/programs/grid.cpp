#include <cstdlib>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <MishMesh/TriMesh.h>
#include <MishMesh/visualization.h>

int main(int argc, char **argv) {
	int resolution[3] = {10, 10, 10};
	if(argc > 1) {
		int res = std::atoi(argv[1]);
		resolution[0] = res;
		resolution[1] = res;
		resolution[2] = res;
	}
	double bbox_ltf[3] = {0,0,0};
	double bbox_rbn[3] = {1,1,1};
	auto mesh = MishMesh::grid_mesh(resolution, bbox_ltf, bbox_rbn, 0.2 / resolution[0]);
	OpenMesh::IO::write_mesh(mesh, "grid.obj");
}
