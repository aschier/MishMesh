#include <iostream>

#include "../thirdparty/ProgramOptions.hxx"

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <MishMesh/TriMesh.h>
#include <MishMesh/smoothing.h>

int main(int argc, char **argv) {
	po::parser parser;
	parser["help"]
		.abbreviation('?')
		.description("print this help screen")
		.callback([&] { std::cout << parser << '\n'; });
	parser["input"]
		.abbreviation('i')
		.description("Input file.")
		.type(po::string);
	parser["output"]
		.abbreviation('o')
		.description("Output file.")
		.type(po::string);
	parser["iterations"]
		.abbreviation('n')
		.description("Number of iterations.")
		.type(po::i32)
		.fallback(1);

	if(argc == 1) {
		std::cerr << parser << std::endl;
	}

	parser.parse(argc, argv);

	if(!parser["input"].available() || !parser["output"].available()) {
		std::cerr << "input (-i) and output (-o) parameters are required." << std::endl;
		exit(1);
	}

	std::string input_filename = parser["input"].get().string;
	std::string output_filename = parser["output"].get().string;
	int iterations = parser["iterations"].get().i32;

	MishMesh::TriMesh mesh;
	OpenMesh::IO::read_mesh(mesh, input_filename);

	MishMesh::smooth_mesh(mesh, parser["iterations"].get().i32);

	OpenMesh::IO::write_mesh(mesh, output_filename);

	return 0;
}