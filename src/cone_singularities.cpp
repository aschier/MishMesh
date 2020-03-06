#include "cone_singularities_impl.h"

#include <MishMesh/laplace.h>
#include <MishMesh/utils.h>
#include <MishMesh/macros.h>

#ifndef NDEBUG
#include <iostream>
#define DEBUG_OUTPUT(x) std::cerr << (x) << std::endl
#else
#define DEBUG_OUTPUT(x)
#endif

/*
 * This is an implementation of the iterative algorithm for finding cone singularities in
 *     Ben-Chen, Mirela, Craig Gotsman, and Guy Bunin. 2008.
 *     "Conformal Flattening by Curvature Prescription and Metric Scaling."
 *     Computer Graphics Forum 27 (2): 449-58. https://doi.org/10.1111/j.1467-8659.2008.01142.x.
 */

using namespace std;

/**
 * Create an indicator vector that is 1 at the index and 0 everywhere else.
 * @param length The length of the vector.
 * @param index the index where the vector should contain a 1.
 * @returns The indicator vector.
 */
inline Eigen::VectorXd indicator_vector(size_t length, size_t index) {
	Eigen::VectorXd indicator(length);
	indicator.setZero();
	indicator[index] = 1.0;
	return indicator;
}

/**
 * Get the initial singularities for a mesh.
 * On meshes with boundary, all boundary vertices are singularities. On other meshes, the vertex with maximum gauss curvature
 * is used when the mesh has a positive euler characteristic, the vertex iwht minimum gauss curvature is used when the mesh has
 * a negative euler characteristic and otherwise an empty set is returned.
 * @param mesh The mesh.
 * @param gauss_curvatures A vector with the gauss curvatures for all mesh vertices.
 * @returns a set of vertices that are the initial singularities.
 */
std::set<size_t> MishMesh::cone_singularities::get_initial_singularities(const MishMesh::TriMesh &mesh, const Eigen::VectorXd &gauss_curvatures) {
	// All boundary vertices are singularities
	set<size_t> singularity_indices{};
	for(const auto &v: mesh.vertices()) {
		if(mesh.is_boundary(v)) {
			singularity_indices.insert(v.idx());
		}
	}
	if(!singularity_indices.empty()) {
		return singularity_indices;
	}

	size_t ec = MishMesh::euler_characteristic(mesh);
	if(ec > 0) {
		size_t curvature_argmax;
		gauss_curvatures.maxCoeff(&curvature_argmax);
		singularity_indices.insert(curvature_argmax);
	} else if(ec < 0) {
		size_t curvature_argmin;
		gauss_curvatures.minCoeff(&curvature_argmin);
		singularity_indices.insert(curvature_argmin);
	}
	return singularity_indices;
}

/**
 * Build the laplace solver for a mesh.
 * @param[out] solver a sparse matrix solver.
 * @param mesh The mesh.
 * @param normalized Set if the normalized laplace should be used.
 * @returns A Eigen::ComputationInfo object for the solver, that can be used to check if there were problems factoring the matrix.
 */
Eigen::ComputationInfo MishMesh::cone_singularities::build_solver(Eigen::SparseLU<Eigen::SparseMatrix<double>> &solver, MishMesh::TriMesh &mesh, bool normalized) {
	Eigen::SparseMatrix<double> laplace = MishMesh::laplace_matrix(mesh, normalized, false);
	solver.compute(laplace);
	Eigen::ComputationInfo info = solver.info();
	if(info != Eigen::ComputationInfo::Success) {
		DEBUG_OUTPUT(std::string("Error factorizing the matrix: ") + solver.lastErrorMessage());
	}
	return info;
}

/**
 * Calculate the curvature flow phi from the original curvature and the new curvature.
 * @param[in] K_orig The original curvature.
 * @param[in] K_new The new curvature.
 * @param[in] laplace_LU A LU factorized mesh laplacian matrix of the mesh.
 * @returns The curvature flow phi.
 */
Eigen::VectorXd MishMesh::cone_singularities::compute_curvature_flow(const Eigen::VectorXd &K_orig, const Eigen::VectorXd &K_new, const Eigen::SparseLU<Eigen::SparseMatrix<double>> &laplace_LU) {
	Eigen::VectorXd phi;
	Eigen::VectorXd rhs = K_new - K_orig;
	phi = laplace_LU.solve(rhs);
	return phi;
}

/**
 * Compute the new curvature of a mesh with cone singularities removed.
 * @param[in] K_orig The original curvature.
 * @param singularity_indices A set of indices of vertices which are cone singularities.
 * @param P_matrix the P matrix from the algorithm. Use the build_P_matrix method to build the matrix.
 * @returns A vector with the new curvature.
 * @see build_P_matrix
 */
Eigen::VectorXd MishMesh::cone_singularities::compute_new_curvature(const Eigen::VectorXd &K_orig, const std::set<size_t> &singularity_indices, const Eigen::SparseMatrix<double> &P_matrix) {
	Eigen::VectorXd K_new(K_orig.rows());
	Eigen::VectorXd K_without_S = K_orig;
	for(auto si: singularity_indices) {
		K_without_S[si] = 0;
	}
	Eigen::SparseLU<Eigen::SparseMatrix<double>> LU2;
	LU2.compute(P_matrix);
	K_new.setZero();
	for(auto si: singularity_indices) {
		Eigen::VectorXd indicator = indicator_vector(K_orig.rows(), si);
		Eigen::VectorXd G_i_hat = LU2.solve(indicator);
		K_new[si] = K_orig[si] + G_i_hat.transpose() * K_without_S;
	}
	return K_new;
}


/**
 * Optimize the curvature of the mesh.
 * @param mesh The mesh.
 * @param K_orig The original curvature of the mesh.
 * @param[out] singularity_indices The indices of the vertices, which are chosen as cone singularites.
 * @param laplace_LU A LU factorized mesh laplacian matrix of the mesh.
 * @param phi_epsilon Convergence criterion for the iteration: When phi_max - phi_min < phi_epsilon, the algorithm is converged.
 * @param max_iterations Stopping criterion: After a maximum of max_iterations, the algorithm is stopped.
 * @param use_optimal_phi If the optimal phi should be used, when the convergence criterion (phi_max - phi_min) gets worse in later iterations.
 * @param mode The ConeAdditionMode, which cones should be added in each iteration.
 * @returns The new curvature, i.e. the K from the last iteration, or the K from the iteration with the best convergence criterion, when use_optimal_phi is set.
 */
Eigen::VectorXd MishMesh::cone_singularities::optimize_curvature(MishMesh::TriMesh &mesh, const Eigen::VectorXd &K_orig, std::set<size_t> &singularity_indices, const Eigen::SparseLU<Eigen::SparseMatrix<double>> &laplace_LU, const double phi_epsilon, const int max_iterations, const bool use_optimal_phi, MishMesh::ConeAdditionMode mode)
{
	double min_phi, max_phi;
	int iter = 0;
	Eigen::VectorXd phi;
	Eigen::VectorXd optimal_K_new;
	double optimal_phi_range = numeric_limits<double>::infinity();
	size_t optimal_num_singularities;
	auto optimal_singularity_indices = singularity_indices;
	Eigen::VectorXd K_new = compute_target_gauss_curvature(K_orig, singularity_indices);
	do {
		phi = compute_curvature_flow(K_orig, K_new, laplace_LU);
		min_phi = phi.minCoeff();
		max_phi = phi.maxCoeff();
		size_t argmin_phi, argmax_phi;
		for(auto si: singularity_indices) {
			phi[si] = numeric_limits<double>::infinity();
		}
		phi.minCoeff(&argmin_phi);
		for(auto si: singularity_indices) {
			phi[si] = -numeric_limits<double>::infinity();
		}
		phi.maxCoeff(&argmax_phi);

		if(mode == ConeAdditionMode::MIN || mode == ConeAdditionMode::BOTH) {
			singularity_indices.insert(argmin_phi);
		}
		if(mode == ConeAdditionMode::MAX || mode == ConeAdditionMode::BOTH) {
			singularity_indices.insert(argmax_phi);
		}

		Eigen::SparseMatrix<double> P_matrix = build_P_matrix(mesh, singularity_indices);
		K_new = compute_new_curvature(K_orig, singularity_indices, P_matrix);
		double phi_range = abs(max_phi - min_phi);
		if(use_optimal_phi) {
			if(phi_range < optimal_phi_range) {
				optimal_phi_range = phi_range;
				optimal_K_new = K_new;
				optimal_num_singularities = singularity_indices.size();
				optimal_singularity_indices = singularity_indices;
			}
		} else {
			optimal_num_singularities = singularity_indices.size();
		}
		DEBUG_OUTPUT(string("current max_phi - min_phi: ") + std::to_string(phi_range));
		if(use_optimal_phi) {
			DEBUG_OUTPUT(string("optimal max_phi - min_phi: ") + std::to_string(optimal_phi_range));
		}
		DEBUG_OUTPUT("");
	} while(max_phi - min_phi > phi_epsilon && ++iter < max_iterations);

	DEBUG_OUTPUT(string("original curvature sum: ") + std::to_string(K_orig.sum()));
	DEBUG_OUTPUT(string("new curvature sum: ") + std::to_string(K_new.sum()));
	DEBUG_OUTPUT(string("difference: ") + std::to_string(K_new.sum() - K_orig.sum()));
	DEBUG_OUTPUT(string("singularities used: ") + std::to_string(optimal_num_singularities));
	if(use_optimal_phi) {
		singularity_indices = optimal_singularity_indices;
		return optimal_K_new;
	} else {
		return K_new;
	}
}

/**
 * Compute the unweightes gauss curvature (2*pi - sum(angles)) for each vertex of a mesh.
 * @param mesh The mesh.
 * @returns A vector with the gauss curvature for each vertex in the mesh.
 */
Eigen::VectorXd MishMesh::cone_singularities::compute_gauss_curvature(const MishMesh::TriMesh &mesh) {
	Eigen::VectorXd result(mesh.n_vertices());
	for(auto vh: mesh.vertices()) {
		double gauss_curvature = 2 * M_PI;
		FOR_CVIH(h_it, vh) {
			gauss_curvature -= mesh.calc_sector_angle(*h_it);
		}
		result[vh.idx()] = gauss_curvature;
	}
	return result;
}

/**
 * Compute a vector with the target gauss curvatures, which 1/(number of cone singularities) at the cone singularities and 0 everywhere else.
 * @param original_curvature The original curvature.
 * @param singularity_indices The indices of the cone singularities.
 * @returns A vector with the target gauss curvature.
 */
Eigen::VectorXd MishMesh::cone_singularities::compute_target_gauss_curvature(const Eigen::VectorXd &original_curvature, const set<size_t> &singularity_indices) {
	Eigen::VectorXd K_new = Eigen::VectorXd::Zero(original_curvature.rows());
	for(auto si: singularity_indices) {
		K_new[si] = original_curvature.sum() / singularity_indices.size();
	}
	return K_new;
}

/**
 * Compute the cone singularities for a mesh.
 * @param mesh The mesh.
 * @param epsilon Convergence criterion for the iteration: When phi_max - phi_min < phi_epsilon, the algorithm is converged.
 * @param max_iterations Stopping criterion: After a maximum of max_iterations, the algorithm is stopped.
 * @param mode The ConeAdditionMode which cones should be added in each iteration.
 * @param initial_singularities A initial set of singularities for the algorithm.
 * @returns A vector with the vertex handles of the cone singularities.
 */
std::vector<MishMesh::TriMesh::VertexHandle> MishMesh::compute_cone_singularities(MishMesh::TriMesh &mesh, double epsilon, int max_iterations, const MishMesh::ConeAdditionMode mode, set<size_t> initial_singularities) {
	Eigen::SparseLU<Eigen::SparseMatrix<double>> laplace_LU;
	MishMesh::cone_singularities::build_solver(laplace_LU, mesh, false);
	Eigen::VectorXd K_orig = MishMesh::cone_singularities::compute_gauss_curvature(mesh);

	set<size_t> singularity_indices;
	if(initial_singularities.size() > 0) {
		// user provided initial singularities
		singularity_indices = initial_singularities;
	} else {
		singularity_indices = MishMesh::cone_singularities::get_initial_singularities(mesh, K_orig);
	}

	Eigen::VectorXd result_curvature = MishMesh::cone_singularities::optimize_curvature(mesh, K_orig, singularity_indices, laplace_LU, epsilon, max_iterations, false, mode);
	std::vector<MishMesh::TriMesh::VertexHandle> result{};
	std::transform(singularity_indices.begin(), singularity_indices.end(), std::back_inserter(result), [&](size_t idx){return mesh.vertex_handle(static_cast<uint>(idx)); });
	return result;
}

/**
 * Build the transition matrix P from the paper.
 * @param mesh The mesh.
 * @param singularity_indices The indices of the cone singularities.
 * @returns The matrix.
 */
Eigen::SparseMatrix<double> MishMesh::cone_singularities::build_P_matrix(MishMesh::TriMesh &mesh, const set<size_t> &singularity_indices) {
	auto triplets = MishMesh::laplace_triplets(mesh, false);
	std::vector<Eigen::Triplet<double>> new_triplets;
	for(int i = 0; i < triplets.size(); i++) {
		if(std::find(singularity_indices.begin(), singularity_indices.end(), triplets[i].row()) != singularity_indices.end()) {
			// triplet belongs to a row in the lower half of the matrix (singularity to ... transitions)
			// so we omit it keeping the lower part zero
			continue;
		} else {
			// Add the upper half of the matrix normally
			new_triplets.push_back(triplets[i]);
		}
	}
	// Fill in the lower right part of the matrix with singularity to singularity transitions with the identity matrix
	for(auto si: singularity_indices) {
		new_triplets.push_back(Eigen::Triplet<double>(static_cast<uint>(si), static_cast<uint>(si), 1.0));
	}
	Eigen::SparseMatrix<double> P(mesh.n_vertices(), mesh.n_vertices());
	P.setFromTriplets(new_triplets.begin(), new_triplets.end());
	return P;
}
