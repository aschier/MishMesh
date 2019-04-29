#pragma once

#include <set>
#include <Eigen/Eigen>
#include <MishMesh/dijkstra.h>

namespace MishMesh {
	std::set<TriMesh::EdgeHandle> minimum_spanning_tree(TriMesh &mesh, std::vector<TriMesh::VertexHandle> vertices, double edge_cost_function(TriMesh &mesh, const TriMesh::HalfedgeHandle edge, const void *param) = edge_length, void *param = nullptr);
}
