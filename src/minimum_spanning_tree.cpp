#include "MishMesh/minimum_spanning_tree.h"
#include <iostream>
#include "MishMesh/dijkstra.h"
#include "MishMesh/search.h"

using namespace std;

namespace MishMesh {
	/**
	 * Compute the minimum spanning trees of the shortest paths between a set of vertices.
	 * @param mesh The mesh.
	 * @param vertices A list of target VertexHandles in the mesh.
	 * @param edge_cost_function A function, that assigns HalfEdgeHandles a cost. By default, the length of the edge is used.
	 * @param edge_cost_param A pointer to additional data passed to edge_cost_function.
	 * @returns A list of edge sets, that form different spanning trees on the mesh connected components of the vertices.
	 */
	template<typename MeshT>
	std::vector<std::set<typename MeshT::EdgeHandle>> minimum_spanning_trees(MeshT &mesh, std::vector<typename MeshT::VertexHandle> vertices, double edge_cost_function(MeshT &mesh, const typename MeshT::HalfedgeHandle edge, const void *param), void *edge_cost_param) {
		vector<vector<typename MeshT::VertexHandle>> vertex_subsets;
		auto cc_vertex_sets = get_connected_components_vertices(mesh);
		for(auto &cc_vertex_set : cc_vertex_sets) {
			vector<typename MeshT::VertexHandle> vertex_subset;
			for(auto vh : vertices) {
				if(cc_vertex_set.find(vh) != cc_vertex_set.end()) {
					vertex_subset.push_back(vh);
				}
			}
			if(!vertex_subset.empty()) {
				vertex_subsets.push_back(vertex_subset);
			}
		}
		vector<set<typename MeshT::EdgeHandle>> result;
		for(auto &vertex_subset : vertex_subsets) {
			const auto mst_result = minimum_spanning_tree(mesh, vertex_subset, edge_cost_function, edge_cost_param);
			result.push_back(mst_result.get_edges());
		}
		return result;
	}

	template<typename MeshT>
	priority_queue<PathEdge<MeshT>, vector<PathEdge<MeshT>>, GreaterPathlength<MeshT>> initialize_search(MeshT &mesh, const OpenMesh::ArrayKernel::VertexHandle &start_vh, double(*edge_cost_function)(MeshT &mesh, OpenMesh::ArrayKernel::HalfedgeHandle edge, const void *param), void *edge_cost_param, OpenMesh::HPropHandleT<double> &prop_edge_shortest_path_length, const OpenMesh::VPropHandleT<double> &prop_vertex_shortest_path_length) {
		// Initialize distance properties
		for(auto &v : mesh.vertices()) {
			mesh.property(prop_vertex_shortest_path_length, v) = numeric_limits<double>::infinity();
		}
		for(auto &h : mesh.halfedges()) {
			mesh.property(prop_edge_shortest_path_length, h) = numeric_limits<double>::infinity();
		}
		mesh.property(prop_vertex_shortest_path_length, start_vh) = 0;

		priority_queue<PathEdge<MeshT>, vector<PathEdge<MeshT>>, GreaterPathlength<MeshT>> queue;
		for(auto h_it = mesh.cvoh_ccwbegin(start_vh); h_it != mesh.cvoh_ccwend(start_vh); h_it++) {
			const auto vh2 = mesh.to_vertex_handle(*h_it);
			double distance = edge_cost_function(mesh, *h_it, edge_cost_param);
			mesh.property(prop_edge_shortest_path_length, *h_it) = distance;
			mesh.property(prop_vertex_shortest_path_length, vh2) = distance;
			queue.push(PathEdge<MeshT>{&mesh, &prop_edge_shortest_path_length, *h_it});
		}
		return queue;
	}

	/**
	 * Calculate all shortest distances from a given vertex on the mesh.
	 * @param mesh The mesh.
	 * @param start_vh The vertex from which the distances are calculated.
	 * @param prop_vertex_shortest_path_length A vertex property, in which the distance from the start vertex is stored.
	 * @param prop_edge_shortest_path_length A half edge property, in which the distance of the to_vertex from the start_vertex is stored.
	 * @param edge_cost_function An edge cost function for measuring path lengths.
	 * @param edge_cost_param A pointer to parameters for the edge cost function.
	 */
	template<typename MeshT>
	void calc_shortest_distances(MeshT &mesh,
		typename MeshT::VertexHandle start_vh,
		OpenMesh::VPropHandleT<double> &prop_vertex_shortest_path_length,
		OpenMesh::HPropHandleT<double> &prop_edge_shortest_path_length,
		double edge_cost_function(MeshT &mesh, const typename MeshT::HalfedgeHandle edge, const void *param),
		void *edge_cost_param) {
		auto queue = initialize_search(mesh, start_vh, edge_cost_function, edge_cost_param, prop_edge_shortest_path_length, prop_vertex_shortest_path_length);
		while(!queue.empty()){
			PathEdge<MeshT> path_edge = queue.top();
			queue.pop();
			auto &heh = path_edge.halfedge_handle;
			auto vh = mesh.to_vertex_handle(heh);
			if(mesh.status(vh).tagged()) continue;
			mesh.status(vh).set_tagged(true);

			// Iterate over all outgoing halfedges of the current vertex
			for(auto h_it = mesh.cvoh_ccwbegin(vh); h_it != mesh.cvoh_ccwend(vh); h_it++) {
				auto to_vh = mesh.to_vertex_handle(*h_it);
				if(mesh.status(to_vh).tagged()) continue;

				double distance = mesh.property(prop_edge_shortest_path_length, heh) + edge_cost_function(mesh, *h_it, edge_cost_param);
				mesh.property(prop_edge_shortest_path_length, *h_it) = distance;
				mesh.property(prop_vertex_shortest_path_length, to_vh) = distance;

				// Add the adjacent unvisited vertices to the queue
				if(!mesh.status(to_vh).tagged()) {
					queue.push(PathEdge<MeshT>{&mesh, &prop_edge_shortest_path_length, *h_it});
				}
			}
		}
	}

	template<typename MeshT>
	MSTResult<MeshT> minimum_spanning_tree(MeshT &mesh, std::vector<typename MeshT::VertexHandle> vertices, double edge_cost_function(MeshT &mesh, const typename MeshT::HalfedgeHandle edge, const void *param), void *edge_cost_param) {
		MSTResult<MeshT> result;
		if(vertices.size() < 2){
			return {};
		}
		set<typename MeshT::VertexHandle> targetVertices(vertices.begin(), vertices.end());
		assert(targetVertices.size() >= 2);

		// Add the needed properties
		OpenMesh::VPropHandleT<double> prop_vertex_shortest_path_length;
		OpenMesh::HPropHandleT<double> prop_edge_shortest_path_length;
		mesh.add_property(prop_vertex_shortest_path_length);
		mesh.add_property(prop_edge_shortest_path_length);
		mesh.request_vertex_status();

		// Use the first vertex as start vertex
		auto start_vh = *targetVertices.begin();
		targetVertices.erase(start_vh);

		// Initialize the queue with the edges reachable from the source vertex
		priority_queue<PathEdge<MeshT>, vector<PathEdge<MeshT>>, GreaterPathlength<MeshT>> queue = initialize_search(mesh, start_vh, edge_cost_function, edge_cost_param, prop_edge_shortest_path_length, prop_vertex_shortest_path_length);

		for(auto &vh : mesh.vertices()) {
			mesh.status(vh).set_tagged(false); // not visited
		}
		mesh.status(start_vh).set_tagged(true);

		/*
		   The edges of the current region store the distance to the closest source vertex.
		   When a new target vertex is found, the path to the closest source vertex is traced backwards.
		   Afterwards the target vertex is added as new source vertex by setting its distance property
		   to 0 and adding all new edges without a previous distance (i.e. their distance is their
		   own length).
		   The already visited edges now do NOT store a correct distance to the nearest source vertex,
		   but because we use a priority queue, all NEW edges will be visited via a path that points
		   to the nearest source vertex. This is true, because the paths from the new source vertex are
		   traced before other paths, until their length becomes longer than the length of previous edges
		   in the priority queue. Then the distance field will continue to grow in all directions.
		   This makes sure, that each newly visited vertex can use the distance field to trace a path to
		   the nearest source vertex, even when the paths between existing source vertices cannot be traced.

		   Explanation: All vertices visited after adding the n-th target vertex to the set of source vertices
		   have the correct distance to the nearest source vertex, because the edges are added using a priority
		   queue. Any vertex affected by a (now) wrong distances was already found and if it is the target vertex,
		   the corresponding path was already added to the result paths. The result is a MST, because if there
		   would be a shorter path from the i-th target vertex (i!=n) to the (n+1)-th target vertex, the (n+1)-th
		   target vertex would have been found before the n-th target vertex.
		 */

		 // As long as there are unfound target vertices search for a new path to any of the target vertices,
		 // until the list of targetVertices is empty.
		while(!targetVertices.empty()) {
			typename MeshT::VertexHandle target_vh;
			bool found_target = false;
			while(!queue.empty() && !found_target) {
				PathEdge<MeshT> path_edge = queue.top();
				queue.pop();
				auto &heh = path_edge.halfedge_handle;
				auto vh = mesh.to_vertex_handle(heh);
				if(mesh.status(vh).tagged()) continue;
				mesh.status(vh).set_tagged(true);

				// Iterate over all outgoing halfedges of the current vertex
				for(auto h_it = mesh.cvoh_ccwbegin(vh); h_it != mesh.cvoh_ccwend(vh); h_it++) {
					auto to_vh = mesh.to_vertex_handle(*h_it);
					if(mesh.status(to_vh).tagged()) continue;

					double distance = mesh.property(prop_edge_shortest_path_length, heh) + edge_cost_function(mesh, *h_it, edge_cost_param);
					mesh.property(prop_edge_shortest_path_length, *h_it) = distance;
					mesh.property(prop_vertex_shortest_path_length, to_vh) = distance;

					// Add the adjacent unvisited vertices to the queue
					if(!mesh.status(to_vh).tagged()) {
						queue.push(PathEdge<MeshT>{&mesh, &prop_edge_shortest_path_length, *h_it});
					}
					// When the vertex is a target vertex, we found a new path in the MST.
					if(targetVertices.find(to_vh) != targetVertices.end()) {
						// Update the distances for the new source vertex and the adjacent edges.
						distance = edge_cost_function(mesh, *h_it, edge_cost_param);
						mesh.property(prop_edge_shortest_path_length, *h_it) = distance;
						mesh.property(prop_vertex_shortest_path_length, to_vh) = 0;
						// Stop the current path search and set the result target_vh
						targetVertices.erase(to_vh);
						target_vh = to_vh;
						found_target = true;
						break;
					}
				}
			};

			if(found_target) {
				// Trace the path backwards from the new target vertex to the closest source vertex.
				DijkstraResult<MeshT> dijkstra_result = trace_path(mesh, target_vh, prop_vertex_shortest_path_length, prop_edge_shortest_path_length);
				result.paths.push_back(dijkstra_result);
			}

			if(!found_target && !targetVertices.empty()) {
				// No target vertex found in this component. Restart search at another targetVertex
				start_vh = *targetVertices.begin();
				targetVertices.erase(targetVertices.begin());
				queue = initialize_search(mesh, start_vh, edge_cost_function, edge_cost_param, prop_edge_shortest_path_length, prop_vertex_shortest_path_length);
				for(auto &vh : mesh.vertices()) {
					mesh.status(vh).set_tagged(false); // not visited
				}
				mesh.status(start_vh).set_tagged(true);
				continue;
			}
		}

		mesh.remove_property(prop_vertex_shortest_path_length);
		mesh.remove_property(prop_edge_shortest_path_length);

		return result;
	}


	template void calc_shortest_distances(TriMesh &mesh,
		typename TriMesh::VertexHandle start_vh,
		OpenMesh::VPropHandleT<double> &prop_vertex_shortest_path_length,
		OpenMesh::HPropHandleT<double> &prop_edge_shortest_path_length,
		double edge_cost_function(TriMesh &mesh, const typename TriMesh::HalfedgeHandle edge, const void *param),
		void *edge_cost_param);
	template void calc_shortest_distances(PolyMesh &mesh,
		typename PolyMesh::VertexHandle start_vh,
		OpenMesh::VPropHandleT<double> &prop_vertex_shortest_path_length,
		OpenMesh::HPropHandleT<double> &prop_edge_shortest_path_length,
		double edge_cost_function(PolyMesh &mesh, const typename PolyMesh::HalfedgeHandle edge, const void *param),
		void *edge_cost_param);

	template MSTResult<TriMesh> minimum_spanning_tree(TriMesh &mesh, std::vector<TriMesh::VertexHandle> vertices, double edge_cost_function(TriMesh &mesh, const TriMesh::HalfedgeHandle edge, const void *param), void *edge_cost_param);
	template MSTResult<PolyMesh> minimum_spanning_tree(PolyMesh &mesh, std::vector<PolyMesh::VertexHandle> vertices, double edge_cost_function(PolyMesh &mesh, const PolyMesh::HalfedgeHandle edge, const void *param), void *edge_cost_param);

	template std::vector<std::set<TriMesh::EdgeHandle>> minimum_spanning_trees(TriMesh &mesh, std::vector<TriMesh::VertexHandle> vertices, double edge_cost_function(TriMesh &mesh, const TriMesh::HalfedgeHandle edge, const void *param), void *edge_cost_param);
	template std::vector<std::set<PolyMesh::EdgeHandle>> minimum_spanning_trees(PolyMesh &mesh, std::vector<PolyMesh::VertexHandle> vertices, double edge_cost_function(PolyMesh &mesh, const PolyMesh::HalfedgeHandle edge, const void *param), void *edge_cost_param);
}
