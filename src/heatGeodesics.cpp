#include "MishMesh/geodesics.h"

#include <MishMesh/laplace.h>
#include <MishMesh/utils.h>
#include <MishMesh/macros.h>

#include <Eigen/Eigen>

#include <vector>
#include <cassert>

/**
 * Compute geodesic distances on a mesh using the heat method by Crane, Wischedel and Wardetzky.
 * [Crane, K., Weischedel, C., and Wardetzky, M. (2013). Geodesics in Heat. ACM Trans. Graph. 32, 1–11.]
 * https://dl.acm.org/citation.cfm?id=2516977
 *
 * @param[inout] mesh The mesh.
 * @param[in] start_vh A valid vertex in the mesh, that will be used as start vertex.
 * @param[in] geodesicGeodesicDistanceProperty A mesh property to store the geodesic distances. The method assumes, that the property is already added to the mesh
 * @param[in] t The timestep for the heat diffusion step in the algorithm.
 */
void MishMesh::compute_heat_geodesics(TriMesh &mesh, const TriMesh::VertexHandle start_vh, GeodesicDistanceProperty geodesicGeodesicDistanceProperty, double t) {
	mesh.request_face_normals();
	mesh.update_face_normals();

	// Lc is the laplacian without area weighting
	Eigen::SparseMatrix<double> Lc = MishMesh::laplace_matrix(mesh, false, false);

	// Set the initial heat distribution
	Eigen::VectorXd u0 = Eigen::VectorXd::Zero(mesh.n_vertices());
	u0[start_vh.idx()] = 1.0;

	// Calculate the vertex areas as 1/3 of the triangles around a vertex
	Eigen::VectorXd vertex_areas(mesh.n_vertices());
	vertex_areas.setZero();
	for(auto vh : mesh.vertices()) {
		FOR_CVF(f_it, vh) {
			vertex_areas[vh.idx()] += MishMesh::compute_area(mesh, *f_it);
		}
		vertex_areas[vh.idx()] /= 3.0;
	}

	// Calculate A - t*L_c
	Eigen::SparseMatrix<double> A_tLc = -t * Lc;
	assert(A_tLc.rows() == A_tLc.cols());
	for(int i = 0; i < A_tLc.rows(); i++){
		A_tLc.coeffRef(i, i) += vertex_areas[i];
	};
	assert(A_tLc.isCompressed());

	// Compute one heat diffusion timestep
	Eigen::SparseLU<Eigen::SparseMatrix<double>> luSolver;
	luSolver.compute(A_tLc);
	Eigen::VectorXd u = luSolver.solve(u0);

	// Calculate the gradients on the faces
	std::vector<OpenMesh::Vec3d> face_grad_u(mesh.n_faces(), OpenMesh::Vec3d(0,0,0));
	for(auto fh : mesh.faces()) {
		OpenMesh::Vec3d &x = face_grad_u[fh.idx()];
		const OpenMesh::Vec3d &N = mesh.normal(fh);
		FOR_CFV(v_it, fh) {
			auto vh = *v_it;
			auto heh = MishMesh::opposite_halfedge(mesh, fh, vh);
			OpenMesh::Vec3d ei = mesh.point(mesh.to_vertex_handle(heh)) - mesh.point(mesh.from_vertex_handle(heh));
			x += u[v_it->idx()] * (N % ei);
		}
		x /= 2 * MishMesh::compute_area(mesh, fh);
	}

	// Calculate the divergence at the vertices
	Eigen::VectorXd vertex_div_u = Eigen::VectorXd::Zero(mesh.n_vertices());
	for(auto vh : mesh.vertices()) {
		auto pi = mesh.point(vh);
		double &div = vertex_div_u[vh.idx()];
		FOR_CVOH(h_it, vh) {
			// The edge and the next edge belong to a common face,
			// except when the current edge is a boundary edge. In that case
			// we skip it, as it will be the next_heh of another halfedge later.
			auto heh = *h_it;
			if(mesh.is_boundary(heh)) continue;
			auto next_h_it = h_it;
			next_h_it++;
			auto next_heh = *next_h_it;

			assert(mesh.face_handle(heh) == mesh.face_handle(mesh.opposite_halfedge_handle(next_heh)));
			auto fh = mesh.face_handle(heh);

			auto p1 = mesh.point(mesh.to_vertex_handle(heh));
			auto p2 = mesh.point(mesh.to_vertex_handle(next_heh));
			OpenMesh::Vec3d e1 = (p1 - pi);
			OpenMesh::Vec3d e2 = (p2 - pi);
			OpenMesh::Vec3d e3 = (p2 - p1);
			OpenMesh::Vec3d X_face = -face_grad_u[fh.idx()] / face_grad_u[fh.idx()].norm();

			double angle1 = acos((-e3.normalized()) | (-e2.normalized()));
			double angle2 = acos(e3.normalized() | (-e1.normalized()));

			div += (1.0 / tan(angle1)) * (e1 | X_face) + (1.0 / tan(angle2)) * (e2 | X_face);
		}
		div /= 2.0;
	}

	// The distance of the source vertex is always 0, so we can use it as boundary condition
	// to get an unique solution of the poisson equation.
	std::vector<std::pair<Eigen::Index, double>> bc{std::make_pair(start_vh.idx(), 0.0)};
	MishMesh::apply_boundary_conditions(Lc, vertex_div_u, bc);
	luSolver.compute(Lc);
	Eigen::VectorXd phi = luSolver.solve(vertex_div_u);
	for(auto vh : mesh.vertices()) {
		mesh.property(geodesicGeodesicDistanceProperty, vh) = phi[vh.idx()];
	}

	mesh.release_face_normals();
}
