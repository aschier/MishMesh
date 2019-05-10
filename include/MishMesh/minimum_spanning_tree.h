#pragma once

#include <set>
#include <Eigen/Eigen>
#include <MishMesh/dijkstra.h>

namespace MishMesh {
	struct MSTResult {
		std::vector<DijkstraResult> paths;
		const std::set<TriMesh::EdgeHandle> get_edges() const {
			std::set<TriMesh::EdgeHandle> result;
			for(auto dr: paths) {
				result.insert(dr.edges.begin(), dr.edges.end());
			}
			return result;
		}
	};

	MSTResult minimum_spanning_tree(TriMesh &mesh, std::vector<TriMesh::VertexHandle> vertices, double edge_cost_function(TriMesh &mesh, const TriMesh::HalfedgeHandle edge, const void *param) = edge_length, void *param = nullptr);
	std::vector<std::set<TriMesh::EdgeHandle>> minimum_spanning_trees(TriMesh &mesh, std::vector<TriMesh::VertexHandle> vertices, double edge_cost_function(TriMesh &mesh, const TriMesh::HalfedgeHandle edge, const void *param), void *edge_cost_param);
}
