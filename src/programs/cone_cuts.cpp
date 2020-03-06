#include <OpenMesh/Core/IO/MeshIO.hh>

#include <MishMesh/TriMesh.h>
#include <MishMesh/cone_singularities.h>
#include <MishMesh/visualization.h>
#include <MishMesh/minimum_spanning_tree.h>

constexpr double LINE_THICKNESS = 0.001;

using namespace std;

int main(int argc, char **argv) {
	MishMesh::TriMesh mesh;
	std::string filename(argv[1]);
	OpenMesh::IO::read_mesh(mesh, filename);

	auto vertices = MishMesh::compute_cone_singularities(mesh, 1.0);
	auto mst = MishMesh::minimum_spanning_tree(mesh, vertices);
	auto cut_edge_set = mst.get_edges();
	vector<MishMesh::TriMesh::EdgeHandle> cut_edges(cut_edge_set.begin(), cut_edge_set.end());

	MishMesh::TriMesh edge_mesh = MishMesh::edge_mesh(mesh, cut_edges, LINE_THICKNESS);

	OpenMesh::IO::write_mesh(edge_mesh, filename + "_cone_tree.obj");

	return 0;
}