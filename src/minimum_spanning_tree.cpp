#include "MishMesh/minimum_spanning_tree.h"
#include <iostream>
#include "MishMesh/dijkstra.h"
#include "MishMesh/search.h"

using namespace std;

namespace MishMesh {
	/**
	 * Compute the minimum spanning tree of the shortest paths between a set of vertices.
	 * @param mesh The mesh.
	 * @param vertices A list of target VertexHandles in the mesh.
	 * @param edge_cost_function A function, that assigns HalfEdgeHandles a cost. By default, the length of the edge is used.
	 * @param edge_cost_param A pointer to additional data passed to edge_cost_function.
	 * @returns A set of edges, that form a minimum spanning tree for the set of vertices.
	 * @note When there are vertices on different connected components of the mesh, only the vertices on one component will be connected.
	 *       Use minimum_spanning_trees to build multiple spanning trees, so that all vertices belong to a spanning tree.
	 */
	std::set<TriMesh::EdgeHandle> minimum_spanning_tree(TriMesh &mesh, std::vector<TriMesh::VertexHandle> vertices, double edge_cost_function(TriMesh &mesh, const TriMesh::HalfedgeHandle edge, const void *param), void *edge_cost_param) {
		set<TriMesh::EdgeHandle> result;

		// Get all paths
		const int num_paths = static_cast<int>((vertices.size() * (vertices.size() - 1)) / 2);
		vector<DijkstraResult> paths(num_paths);
		int idx = 0;
#pragma omp parallel for
		for(int i = 0; i < vertices.size(); i++) {
#ifdef _OPENMP
			// copy mesh when running dijkstra in parallel, because the temporarily changes mesh properties
			auto mymesh = mesh;
#else
			auto &mymesh = mesh;
#endif
			for(int j = i + 1; j < vertices.size(); j++) {
				auto result = dijkstra(vertices[i], vertices[j], mymesh, edge_cost_function, edge_cost_param);
				if(result.is_valid()) {
					paths[idx] = result;
					idx++;
				}
			}
		}

		set<TriMesh::VertexHandle> tree_vertices{vertices[0]};
		uint count = 0;
		// In each step add the shortest path, that connects a new target vertex to the graph
		while(tree_vertices.size() < vertices.size() && count <= paths.size()) {
			DijkstraResult best_path;
			best_path.length = numeric_limits<double>::infinity();
			for(auto path: paths) {
				if(!path.is_valid()) continue;
				assert(path.vertices.size() > 0);
				bool v1_in_graph = tree_vertices.find(path.vertices.front()) != tree_vertices.end();
				bool v2_in_graph = tree_vertices.find(path.vertices.back()) != tree_vertices.end();
				if(v1_in_graph != v2_in_graph && path.length < best_path.length) {
					best_path = path;
				}
			}
			if(best_path.length < numeric_limits<double>::infinity()) {
				result.insert(best_path.edges.begin(), best_path.edges.end());
				// either the start or the target vertex is already in the set
				tree_vertices.insert(best_path.vertices.front());
				tree_vertices.insert(best_path.vertices.back());
			}
			count++;
		}
		// The loop terminated before all vertices were added to the tree, because the maximum number of possible paths was tested.
		if(count > paths.size()) {
			std::cerr << "Some vertices could not be added to the spanning tree, because there is no path with finite cost." << endl;
		}
		return result;
	}
	/**
	 * Compute the minimum spanning trees of the shortest paths between a set of vertices.
	 * @param mesh The mesh.
	 * @param vertices A list of target VertexHandles in the mesh.
	 * @param edge_cost_function A function, that assigns HalfEdgeHandles a cost. By default, the length of the edge is used.
	 * @param edge_cost_param A pointer to additional data passed to edge_cost_function.
	 * @returns A list of edge sets, that form different spanning trees on the mesh connected components of the vertices.
	 */
	std::vector<std::set<TriMesh::EdgeHandle>> minimum_spanning_trees(TriMesh &mesh, std::vector<TriMesh::VertexHandle> vertices, double edge_cost_function(TriMesh &mesh, const TriMesh::HalfedgeHandle edge, const void *param), void *edge_cost_param) {
		vector<vector<TriMesh::VertexHandle>> vertex_subsets;
		auto cc_vertex_sets = get_connected_components_vertices(mesh);
		for(auto &cc_vertex_set : cc_vertex_sets) {
			vector<TriMesh::VertexHandle> vertex_subset;
			for(auto vh : vertices) {
				if(cc_vertex_set.find(vh) != cc_vertex_set.end()) {
					vertex_subset.push_back(vh);
				}
			}
			if(!vertex_subset.empty()) {
				vertex_subsets.push_back(vertex_subset);
			}
		}
		vector<set<TriMesh::EdgeHandle>> result;
		for(auto &vertex_subset : vertex_subsets) {
			const auto edges = minimum_spanning_tree(mesh, vertex_subset, edge_cost_function, edge_cost_param);
			result.push_back(edges);
		}
		return result;
	}
}
