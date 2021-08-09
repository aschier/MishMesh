#include <OpenMesh/Core/IO/MeshIO.hh>
#include <MishMesh/TriMesh.h>
#include <MishMesh/utils.h>

int main(int argc, char **argv) {
	if(argc == 1) exit(1);

	MishMesh::TriMesh mesh;
	OpenMesh::IO::read_mesh(mesh, argv[1]);
	MishMesh::Flatness flatness = MishMesh::is_flat(mesh);

	switch(flatness) {
	case MishMesh::Flatness::X:
		std::cout << "Flat: Constant X coordinates." << std::endl;
		break;
	case MishMesh::Flatness::Y:
		std::cout << "Flat: Constant Y coordinates." << std::endl;
		break;
	case MishMesh::Flatness::Z:
		std::cout << "Flat: Constant Z coordinates." << std::endl;
		break;
	case MishMesh::Flatness::NONFLAT:
		std::cout << "Not flat." << std::endl;
		break;
	case MishMesh::Flatness::DEGENERATE:
		std::cout << "Degenerate mesh." << std::endl;
		break;
	}
}