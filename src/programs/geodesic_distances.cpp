#include <iostream>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <MishMesh/TriMesh.h>
#include <MishMesh/visualization.h>
#include <MishMesh/geodesics.h>

#include "../thirdparty/ProgramOptions.hxx"

int main(int argc, char **argv) {
	po::parser parser;
	parser["help"]
		.abbreviation('h')
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
	parser["startvertex"]
		.abbreviation('s')
		.type(po::i32)
		.description("The start vertex for calculating geodesic distances.")
		.fallback(0);
#ifdef HAS_EIGEN
	parser["method"]
		.abbreviation('m')
		.type(po::string)
		.description("The method to compute geodesics distances. Options: \"novotni\" (\"n\"), \"heat\" (\"h\").")
		.fallback("n");
	parser["timestep"]
		.abbreviation('t')
		.type(po::f32)
		.description("The timestep used for the Geodesics in Heat method.")
		.fallback(0.1);
#endif

	parser.parse(argc, argv);

	if(!parser["input"].available() || !parser["output"].available()) {
		std::cerr << "Please specify an input file and an output filename prefix." << std::endl;
		return 1;
	}

	MishMesh::TriMesh mesh;
	OpenMesh::IO::read_mesh(mesh, parser["input"].get().string);

	MishMesh::GeodesicDistanceProperty distanceProperty;
	mesh.add_property(distanceProperty);

#ifdef HAS_EIGEN
	if(parser["method"].get().string == "novotni" || parser["method"].get().string == "n") {
		MishMesh::compute_novotni_geodesics(mesh, mesh.vertex_handle(parser["startvertex"].get().i32), distanceProperty);
	} else if(parser["method"].get().string == "heat" || parser["method"].get().string == "h") {
		MishMesh::compute_heat_geodesics(mesh, mesh.vertex_handle(parser["startvertex"].get().i32), distanceProperty, parser["timestep"].get().f32);
	} else {
		std::cerr << "Unknown method." << std::endl;
		exit(1);
	}
#else
	MishMesh::compute_novotni_geodesics(mesh, mesh.vertex_handle(parser["startvertex"].get().i32), distanceProperty);
#endif

	mesh.request_vertex_colors();
	MishMesh::colorize_mesh(mesh, distanceProperty);
	OpenMesh::IO::write_mesh(mesh, parser["output"].get().string, OpenMesh::IO::Options::VertexColor);

	return 0;
}
