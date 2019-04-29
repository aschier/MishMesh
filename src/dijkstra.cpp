#include "MishMesh/dijkstra.h"
#include <iostream>

using namespace std;

namespace MishMesh {
	/**
	 * An edge on a path. PathEdges can be ordered by the length of the path up to this edge.
	 * @param mesh The mesh the edge belongs to. Used for length calculations.
	 * @param The property handle for the "shortest path" property on the mesh.
	 * @param halfedge_handle The handle to the edge in the mesh.
	 */
	class PathEdge {
	public:
		TriMesh *mesh;
		OpenMesh::HPropHandleT<double> *prop_edge_shortest_path_length;
		TriMesh::HalfedgeHandle halfedge_handle;
	};
	/**
	 * A comparator, that compares two PathEdges by the length of the path up to this edge.
	 */
	class GreaterPathlengh {
	public:
		bool operator()(const PathEdge &a, const PathEdge &b) const {
			return a.mesh->property(*a.prop_edge_shortest_path_length, a.halfedge_handle) > b.mesh->property(*b.prop_edge_shortest_path_length, b.halfedge_handle);
		}
	};

	/**
	 * Find the shortest path from a start vertex to a target vertex using Dijkstra's algorithm.
	 * @param start_vh A VertexHandle for the start vertex.
	 * @param target_vh A VertexHandle for the target vertex.
	 * @param mesh The mesh.
	 * @param edge_cost_function A function, that assigns a cost to a HalfedgeHandle. By default, the length of the edge is used.
	 * @param edge_cost_param A pointer to additional data passed to edge_cost_function.
	 */
	DijkstraResult dijkstra(const TriMesh::VertexHandle start_vh, const TriMesh::VertexHandle target_vh, TriMesh &mesh, double edge_cost_function(TriMesh &mesh, const TriMesh::HalfedgeHandle edge, const void *param), void *edge_cost_param) {
		DijkstraResult result;

		OpenMesh::VPropHandleT<double> prop_vertex_shortest_path_length;
		OpenMesh::HPropHandleT<double> prop_edge_shortest_path_length;
		mesh.add_property(prop_vertex_shortest_path_length);
		mesh.add_property(prop_edge_shortest_path_length);

		for(auto &v : mesh.vertices()) {
			mesh.property(prop_vertex_shortest_path_length, v) = numeric_limits<double>::infinity();
		}
		for(auto &h : mesh.halfedges()) {
			mesh.property(prop_edge_shortest_path_length, h) = numeric_limits<double>::infinity();
		}

		// Initialize the queue with the edges reachable from the source vertex
		priority_queue<PathEdge, vector<PathEdge>, GreaterPathlengh> queue;
		set<TriMesh::VertexHandle> visited_vertices{start_vh};
		mesh.property(prop_vertex_shortest_path_length, start_vh) = 0;
		for(auto h_it = mesh.cvoh_ccwbegin(start_vh); h_it != mesh.cvoh_ccwend(start_vh); h_it++) {
			mesh.property(prop_edge_shortest_path_length, *h_it) = edge_cost_function(mesh, *h_it, edge_cost_param);
			queue.push(PathEdge{&mesh, &prop_edge_shortest_path_length, *h_it});
		}

		do {
			PathEdge path_edge = queue.top();
			queue.pop();
			auto &heh = path_edge.halfedge_handle;
			auto vh = mesh.to_vertex_handle(heh);
			if(visited_vertices.find(vh) != visited_vertices.end()) continue;
			visited_vertices.insert(vh);

			mesh.property(prop_vertex_shortest_path_length, vh) = mesh.property(prop_edge_shortest_path_length, heh) + edge_cost_function(mesh, heh, edge_cost_param);
			// Iterate over all outgoin halfedges of the current vertex
			for(auto h_it = mesh.cvoh_ccwbegin(vh); h_it != mesh.cvoh_ccwend(vh); h_it++) {
				auto vh2 = mesh.to_vertex_handle(*h_it);
				if(visited_vertices.find(vh2) != visited_vertices.end()) continue;
				mesh.property(prop_edge_shortest_path_length, *h_it) = mesh.property(prop_edge_shortest_path_length, heh) + edge_cost_function(mesh, *h_it, edge_cost_param);
				if(vh2 == target_vh) {
					// Empty the queue to stop the outer loop, then break the inner loop.
					queue = {};
					break;
				}
				if(visited_vertices.find(vh2) == visited_vertices.end()) {
					queue.push(PathEdge{&mesh, &prop_edge_shortest_path_length, *h_it});
				}
			}
		} while(!queue.empty());

		// Trace the path backwards
		auto vh = target_vh;
		result.length = mesh.property(prop_vertex_shortest_path_length, target_vh);
		if(result.length == numeric_limits<double>::infinity()) {
			// No path found
			return {};
		}
		result.vertices.push_back(target_vh);
		size_t edge_count = 0;
		size_t max_edge_count = mesh.n_edges();
		do {
			double smallest_distance = numeric_limits<double>::infinity();
			TriMesh::HalfedgeHandle heh;
			for(auto h_it = mesh.cvih_ccwbegin(vh); h_it != mesh.cvih_ccwend(vh); h_it++) {
				double distance = mesh.property(prop_edge_shortest_path_length, *h_it);
				if(distance < smallest_distance) {
					smallest_distance = distance;
					heh = *h_it;
				}
			}
			assert(smallest_distance != numeric_limits<double>::infinity());
			assert(heh.is_valid());
			if(!heh.is_valid()) return{};

			result.edges.push_back(mesh.edge_handle(heh));
			result.vertices.push_back(mesh.from_vertex_handle(heh));
			vh = mesh.from_vertex_handle(heh);

			if(++edge_count > max_edge_count) {
				return {};
			}
		} while(vh != start_vh);
		std::reverse(result.edges.begin(), result.edges.end());
		std::reverse(result.vertices.begin(), result.vertices.end());

		mesh.remove_property(prop_vertex_shortest_path_length);
		mesh.remove_property(prop_edge_shortest_path_length);

		return result;
	}
}
