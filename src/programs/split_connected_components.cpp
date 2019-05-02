#include <iostream>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <MishMesh/TriMesh.h>
#include <MishMesh/split.h>

#include "../thirdparty/ProgramOptions.hxx"

int main(int argc, char **argv) {
	po::parser parser;
	parser["help"]
		.abbreviation('?')
		.description("print this help screen")
		.callback([&] { std::cout << parser << '\n'; });
	parser["input"]
		.abbreviation('i')
		.type(po::string)
		.description("The input file.");
	parser["output"]
		.abbreviation('o')
		.type(po::string)
		.description("The output file prefix.");
	parser.parse(argc, argv);
	if(!parser["input"].available() || !parser["output"].available()) {
		std::cerr << "Please specify an input file and an output filename prefix." << std::endl;
		return 1;
	}

	MishMesh::TriMesh input_mesh;
	OpenMesh::IO::read_mesh(input_mesh, parser["input"].get().string);
	std::vector<MishMesh::TriMesh> output_meshes = MishMesh::split_connected_components(input_mesh);
	int count = 0;
	for(auto &m : output_meshes) {
		OpenMesh::IO::write_mesh(m, parser["output"].get().string + std::to_string(++count) + ".obj");
	}
	return 0;
}
