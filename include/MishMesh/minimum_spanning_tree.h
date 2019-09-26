#pragma once

#include <set>
#include <MishMesh/dijkstra.h>
#include <MishMesh/PolyMesh.h>

namespace MishMesh {
	template<typename MeshT>
	struct MSTResult {
		std::vector<DijkstraResult<MeshT>> paths;
		const std::set<PolyMesh::EdgeHandle> get_edges() const {
			std::set<PolyMesh::EdgeHandle> result;
			for(auto dr: paths) {
				result.insert(dr.edges.begin(), dr.edges.end());
			}
			return result;
		}
	};

	template<typename MeshT>
	MSTResult<MeshT> minimum_spanning_tree(MeshT &mesh, std::vector<typename MeshT::VertexHandle> vertices, double edge_cost_function(MeshT &mesh, const typename MeshT::HalfedgeHandle edge, const void *param) = edge_length, void *param = nullptr);

	template<typename MeshT>
	std::vector<std::set<typename MeshT::EdgeHandle>> minimum_spanning_trees(MeshT &mesh, std::vector<typename MeshT::VertexHandle> vertices, double edge_cost_function(MeshT &mesh, const typename MeshT::HalfedgeHandle edge, const void *param), void *edge_cost_param);
}
