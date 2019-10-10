#include "MishMesh/laplace.h"

#include "MishMesh/macros.h"
#include "MishMesh/utils.h"

#include <cmath>

using namespace std;
using namespace MishMesh;

inline double cot(double angle) {
	return std::cos(angle) / std::sin(angle);
}

/**
 * Create Eigen::Triplet sparse matrix entries for a Laplace matrix.
 * @param mesh The mesh for which the laplacian is built.
 * @param normalized Set if the Laplace matrix should be normalized to have -1 on the diagonal.
 * @param area_weighted When area_weighted is set, the entries are divided by the vertex area (1/3 of the area of the adjacent faces).
 * @note You usually should use MishMesh::laplace_matrix instead of this function, except when you
 *       need to build more complicated matrices, e.g., a block matrix with a laplacian in it.
 */
std::vector<Eigen::Triplet<double>> MishMesh::laplace_triplets(TriMesh &mesh, bool normalized, const bool area_weighted) {
	// the first n_vertices entries are the diagonal entries
	vector<Eigen::Triplet<double>> triplets;
	for(int i = 0; i < mesh.n_vertices(); i++) {
		triplets.push_back(Eigen::Triplet<double>(i, i, 0.0));
	}
	OpenMesh::HPropHandleT<double> prop_halfedge_angle;
	mesh.add_property(prop_halfedge_angle);

	for(auto heh : mesh.halfedges()) {
		if(mesh.is_boundary(heh)) {
			continue;
		}
		auto next_heh = mesh.next_halfedge_handle(heh);
		mesh.property(prop_halfedge_angle, heh) = mesh.calc_sector_angle(next_heh);
		assert(mesh.property(prop_halfedge_angle, heh) <= M_PI);
	}

	for(auto vh : mesh.vertices()) {
		double vertex_area = 0.0;
		if(area_weighted) {
			FOR_CVF(f_it, vh) {
				vertex_area += compute_area(mesh, *f_it);
			}
			vertex_area /= 3.0;
		}
		FOR_CVOH(h_it, vh) {
			const double length = mesh.calc_edge_length(*h_it);
			auto heh1 = *h_it;
			auto heh2 = mesh.opposite_halfedge_handle(*h_it);
			int from_idx = mesh.from_vertex_handle(heh1).idx();
			int to_idx = mesh.to_vertex_handle(heh1).idx();
			double weight = 0.0;
			if(!mesh.is_boundary(heh1)) {
				weight += cot(mesh.property(prop_halfedge_angle, heh1));
			}
			if(!mesh.is_boundary(heh2)) {
				weight += cot(mesh.property(prop_halfedge_angle, heh2));
			}
			weight /= 2.0;
			if(area_weighted) {
				weight /= vertex_area;
			}
			assert(weight != 0);
			triplets[from_idx] = Eigen::Triplet<double>(from_idx, from_idx, triplets[from_idx].value() - weight);
			triplets.push_back(Eigen::Triplet<double>(from_idx, to_idx, weight));
		}
	}

	if(normalized) {
		vector<double> diagonal;
		std::transform(triplets.begin(), triplets.begin() + mesh.n_vertices(), std::back_inserter(diagonal), [&](Eigen::Triplet<double> &t) -> auto {return abs(t.value()); });
		// normalized laplacian = D^{-0.5} L D^{-0.5}
		for(int i = 0; i < triplets.size(); i++) {
			triplets[i] = Eigen::Triplet<double>(triplets[i].row(), triplets[i].col(), triplets[i].value() / sqrt(diagonal[triplets[i].row()] * diagonal[triplets[i].col()]));
		}
	}

	mesh.remove_property(prop_halfedge_angle);
	return triplets;
}

/**
 * Create a sparse Laplace matrix.
 * @param mesh The mesh for which the laplacian is built.
 * @param normalized Set if the Laplace matrix should be normalized to have -1 on the diagonal.
 * @param area_weighted When area_weighted is set, the entries are divided by the vertex area (1/3 of the area of the adjacent faces).
 */
Eigen::SparseMatrix<double> MishMesh::laplace_matrix(TriMesh & mesh, bool normalized, const bool area_weighted) {
	auto triplets = laplace_triplets(mesh, normalized, area_weighted);
	Eigen::SparseMatrix<double> laplace(mesh.n_vertices(), mesh.n_vertices());
	laplace.setFromTriplets(triplets.begin(), triplets.end());
	return laplace;
}

/**
 * Apply dirichlet boundary conditions to a matrix by adding the result for the boundary values to the
 * right hand side and replacing rows/columns of the boundary values in the laplace matrix with 0 .. 0 1 0 .. 0.
 * @param[inout] laplacian The laplace matrix.
 * @param[inout] the right hand side.
 * @param[in] boundary_conditions (index, value) pairs for the boundary conditions.
 * @note The sparsity structure of the matrix is not changed, i.e. the non-zero off-diagonal entries that are
 *       set to zero are explicitly replaced with zero. You need to use makeCompressed on the result matrix,
 *       when you want to compress the zero values.
 */
void MishMesh::apply_boundary_conditions(Eigen::SparseMatrix<double> &laplacian, Eigen::VectorXd &rhs, const vector<pair<Eigen::Index, double>> &boundary_conditions) {
	set<Eigen::Index> boundary_indices;
	for(auto bc : boundary_conditions) {
		boundary_indices.insert(bc.first);
		rhs -= laplacian.col(bc.first) * bc.second;
	}
	for(auto bc : boundary_conditions) {
		rhs[bc.first] = bc.second;
	}
	for(int k = 0; k < laplacian.outerSize(); ++k) {
		for(Eigen::SparseMatrix<double>::InnerIterator it(laplacian, k); it; ++it) {
			if(boundary_indices.find(it.col()) != boundary_indices.end() || boundary_indices.find(it.row()) != boundary_indices.end()) {
				if(it.row() == it.col()) {
					it.valueRef() = 1.0;
				} else {
					it.valueRef() = 0.0;
				}
			}
		}
	}
}
