#pragma once

#include <Eigen/Eigen>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <set>

#include <MishMesh/TriMesh.h>

namespace MishMesh {
	enum ConeAdditionMode {
		MIN,
		MAX,
		BOTH
	};
	std::vector<MishMesh::TriMesh::VertexHandle> compute_cone_singularities(MishMesh::TriMesh &mesh, double epsilon = 1.0, int max_iterations = 10, const ConeAdditionMode mode = ConeAdditionMode::BOTH, std::set<size_t> initial_singularities = {});
}
