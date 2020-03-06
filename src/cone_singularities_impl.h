#pragma once
#include <MishMesh/cone_singularities.h>
#include <MishMesh/TriMesh.h>

#include <Eigen/Eigen>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <set>


namespace MishMesh {
	namespace cone_singularities {
		std::set<size_t> get_initial_singularities(const MishMesh::TriMesh &mesh, const Eigen::VectorXd &K_orig);
		Eigen::ComputationInfo build_solver(Eigen::SparseLU<Eigen::SparseMatrix<double>> &solver, MishMesh::TriMesh &mesh, const bool normalized);

		Eigen::VectorXd compute_curvature_flow(const Eigen::VectorXd &K_orig, const Eigen::VectorXd &K_new, const Eigen::SparseLU<Eigen::SparseMatrix<double>> &laplace_LU);
		Eigen::VectorXd optimize_curvature(MishMesh::TriMesh &mesh, const Eigen::VectorXd &K_orig, std::set<size_t> &singularity_indices, const Eigen::SparseLU<Eigen::SparseMatrix<double>> &laplace_LU, const double phi_epsilon = 1.0, const int max_iterations = 10, const bool use_optimal_phi = false, const ConeAdditionMode mode = ConeAdditionMode::BOTH);
		Eigen::VectorXd compute_new_curvature(const Eigen::VectorXd &K_orig, const std::set<size_t> &singularity_indices, const Eigen::SparseMatrix<double> &P_matrix);
		Eigen::VectorXd compute_gauss_curvature(const MishMesh::TriMesh &mesh);
		Eigen::VectorXd compute_target_gauss_curvature(const Eigen::VectorXd &K_orig, const std::set<size_t> &singularity_indices);
		Eigen::SparseMatrix<double> build_P_matrix(MishMesh::TriMesh &mesh, const std::set<size_t> &singularity_indices);
	}
}