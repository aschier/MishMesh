#pragma once

#include <set>
#include <vector>

#include <MishMesh/TriMesh.h>

#include <Eigen/Eigen>
#include <Eigen/SparseCore>

namespace MishMesh {
	std::vector<Eigen::Triplet<double>> cotan_laplace_triplets(TriMesh &mesh, const bool area_weighted = false);
	void normalize_laplace_triplets(std::vector<Eigen::Triplet<double>> &triplets, const uint rows);
	Eigen::SparseMatrix<double> laplace_matrix(TriMesh &mesh, bool normalized, const bool area_weighted = false);

	struct BoundaryCondition {
		Eigen::Index index;
		double value;
		double diagonal_value = 1.0;
	};
	void apply_boundary_conditions(Eigen::SparseMatrix<double> &laplacian, Eigen::VectorXd &rhs, const std::vector<BoundaryCondition> &boundary_conditions);
	void apply_boundary_conditions(Eigen::SparseMatrix<double> &laplacian, Eigen::VectorXd &rhs, const std::vector<std::pair<Eigen::Index, double>> &boundary_conditions);
}