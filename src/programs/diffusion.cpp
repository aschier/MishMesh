#include <iostream>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <MishMesh/TriMesh.h>
#include <MishMesh/visualization.h>
#include <MishMesh/laplace.h>

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
	parser["steps"]
		.abbreviation('S')
		.type(po::i32)
		.description("The number of diffusion steps.")
		.fallback(1);
	parser["diffusion"]
		.abbreviation('d')
		.type(po::f32)
		.description("The diffusion factor d. The program solves for (I - d*L)*x = rhs.")
		.fallback(10.0);

	parser.parse(argc, argv);

	if(!parser["input"].available() || !parser["output"].available()) {
		std::cerr << "Please specify an input file and an output filename prefix." << std::endl;
		return 1;
	}

	MishMesh::TriMesh mesh;
	OpenMesh::IO::read_mesh(mesh, parser["input"].get().string);

	auto laplace = MishMesh::laplace_matrix(mesh, false, true);
	Eigen::SparseMatrix<double> eye(mesh.n_vertices(), mesh.n_vertices());
	eye.setIdentity();
	Eigen::VectorXd rhs(mesh.n_vertices());
	rhs.setZero();
	int idx = parser["startvertex"].get().i32;
	double value = 1.0;
	rhs[idx] = value;
	double d = parser["diffusion"].get().f32;
	Eigen::SparseMatrix<double> A = eye - d*laplace;

	Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
	solver.compute(A);
	int steps = parser["steps"].get().i32;
	Eigen::VectorXd result;
	for(int i = 0; i < steps; i++) {
		result = solver.solve(rhs);
		rhs = result;
	}
	OpenMesh::VPropHandleT<double> heatProperty;
	mesh.add_property(heatProperty);
	for(auto vh : mesh.vertices()) {
		std::cerr << result[vh.idx()] << std::endl;
		mesh.property(heatProperty, vh) = result[vh.idx()];
	}

	mesh.request_vertex_colors();
	MishMesh::colorize_mesh(mesh, heatProperty);
	OpenMesh::IO::write_mesh(mesh, parser["output"].get().string, OpenMesh::IO::Options::VertexColor);

	return 0;
}
