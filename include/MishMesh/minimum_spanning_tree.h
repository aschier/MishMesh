#pragma once

#include <set>
#include <MishMesh/dijkstra.h>

namespace MishMesh {
	template<typename MeshT>
	struct MSTResult {
		std::vector<DijkstraResult<MeshT>> paths;
		const std::set<typename MeshT::EdgeHandle> get_edges() const {
			std::set<typename MeshT::EdgeHandle> result;
			for(auto dr: paths) {
				result.insert(dr.edges.begin(), dr.edges.end());
			}
			return result;
		}
	};

	template<typename MeshT>
	void calc_shortest_distances(MeshT & mesh, typename MeshT::VertexHandle vertex, OpenMesh::VPropHandleT<double>& prop_vertex_shortest_path_length, OpenMesh::HPropHandleT<double>& prop_edge_shortest_path_length, double edge_cost_function(MeshT &mesh, const typename MeshT::HalfedgeHandle edge, const void *param) = edge_length, void *edge_cost_param = nullptr);

	template<typename MeshT>
	MSTResult<MeshT> minimum_spanning_tree(MeshT &mesh, std::vector<typename MeshT::VertexHandle> vertices, double edge_cost_function(MeshT &mesh, const typename MeshT::HalfedgeHandle edge, const void *param) = edge_length, void *param = nullptr);

	template<typename MeshT>
	std::vector<std::set<typename MeshT::EdgeHandle>> minimum_spanning_trees(MeshT &mesh, std::vector<typename MeshT::VertexHandle> vertices, double edge_cost_function(MeshT &mesh, const typename MeshT::HalfedgeHandle edge, const void *param), void *edge_cost_param);
}
