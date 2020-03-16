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

	namespace cone_singularities {
		typedef Eigen::SparseLU<Eigen::SparseMatrix<double>> LaplaceSolver;
		std::set<size_t> get_initial_singularities(const MishMesh::TriMesh &mesh, const Eigen::VectorXd &K_orig);
		Eigen::ComputationInfo build_solver(Eigen::SparseLU<Eigen::SparseMatrix<double>> &solver, MishMesh::TriMesh &mesh, const bool normalized = false);

		std::vector<MishMesh::TriMesh::VertexHandle> compute_cone_singularities(MishMesh::TriMesh &mesh, LaplaceSolver &solver, double epsilon = 1.0, int max_iterations = 10, const ConeAdditionMode mode = ConeAdditionMode::BOTH, std::set<size_t> initial_singularities = {});
	}
}
