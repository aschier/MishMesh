#pragma once

#include <vector>
#include <set>

#include <MishMesh/TriMesh.h>

#include <Eigen/Eigen>
#include <Eigen/SparseCore>

namespace MishMesh{
	std::vector<Eigen::Triplet<double>> laplace_triplets(TriMesh &mesh, bool normalized, const bool area_weighted = false);
	Eigen::SparseMatrix<double> laplace_matrix(TriMesh &mesh, bool normalized, const bool area_weighted = false);

	struct BoundaryCondition {
		Eigen::Index index;
		double value;
		double diagonal_value = 1.0;
	};
	void apply_boundary_conditions(Eigen::SparseMatrix<double> &laplacian, Eigen::VectorXd &rhs, const std::vector<BoundaryCondition> &boundary_conditions);
	void apply_boundary_conditions(Eigen::SparseMatrix<double> &laplacian, Eigen::VectorXd &rhs, const std::vector<std::pair<Eigen::Index, double>> &boundary_conditions);
}