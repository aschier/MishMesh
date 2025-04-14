#include "MishMesh/utils.h"
#include <MishMesh/macros.h>

#include <MishMesh/PolyMesh.h>
#include <MishMesh/TriMesh.h>

using namespace std;

namespace MishMesh {
	/**
	 * Get the halfedge of a triangle, that is opposite to a given vertex.
	 * @param mesh The mesh.
	 * @param fh The face.
	 * @param vh The vertex.
	 * @returns The halfedge of the face, that is opposite to the vertex.
	 */
	TriMesh::HalfedgeHandle opposite_halfedge(const TriMesh &mesh, const TriMesh::FaceHandle &fh, const TriMesh::VertexHandle &vh) {
		TriMesh::HalfedgeHandle result_heh;
		FOR_CFH(h_it, fh) {
			if(mesh.from_vertex_handle(*h_it) != vh && mesh.to_vertex_handle(*h_it) != vh) {
				result_heh = *h_it;
				break;
			}
		}
		assert(result_heh.is_valid());
		return result_heh;
	}

	/**
	 * Get the vertex handles of a half edge
	 * @param mesh The mesh.
	 * @param heh The halfedge.
	 * @returns The vertex handles of the vertices of the given halfedge.
	 */
	template<typename MeshT>
	std::array<typename MeshT::VertexHandle, 2> edge_vertices(const MeshT &mesh, const typename MeshT::HalfedgeHandle heh) {
		return {mesh.from_vertex_handle(heh), mesh.to_vertex_handle(heh)};
	}

	/**
	 * Get the vertex handles of an edge
	 * @param mesh The mesh.
	 * @param eh The edge.
	 * @returns The vertex handles of the vertices of the given edge.
	 */
	template<typename MeshT>
	std::array<typename MeshT::VertexHandle, 2> edge_vertices(const MeshT &mesh, const typename MeshT::EdgeHandle eh) {
		return edge_vertices(mesh, mesh.halfedge_handle(eh, 0));
	}

	/**
	 * Get the points of a half edge
	 * @param mesh The mesh.
	 * @param heh The halfedge.
	 * @returns The points of the vertices of the given halfedge.
	 */
	template<typename MeshT>
	std::array<typename MeshT::Point, 2> edge_points(const MeshT &mesh, const typename MeshT::HalfedgeHandle heh) {
		return {mesh.point(mesh.from_vertex_handle(heh)), mesh.point(mesh.to_vertex_handle(heh))};
	}

	/**
	 * Get the points of an edge
	 * @param mesh The mesh.
	 * @param eh The edge.
	 * @returns The points of the vertices of the given halfedge.
	 */
	template<typename MeshT>
	std::array<typename MeshT::Point, 2> edge_points(const MeshT &mesh, const typename MeshT::EdgeHandle eh) {
		return edge_points(mesh, mesh.halfedge_handle(eh, 0));
	}

	/**
	 * Get the vertex handles of a face.
	 * @param mesh The mesh.
	 * @param fh The face.
	 * @returns The vertex handles of the vertices of the given face.
	 */
	std::array<TriMesh::VertexHandle, 3> face_vertices(const TriMesh &mesh, const TriMesh::FaceHandle fh) {
		array<TriMesh::VertexHandle, 3> vhs;
		std::copy(mesh.cfv_ccwbegin(fh), mesh.cfv_ccwend(fh), vhs.begin());
		return vhs;
	}

	/**
	 * Get the vertex handles of a face.
	 * @param mesh The mesh.
	 * @param fh The face.
	 * @returns The vertex handles of the vertices of the given face.
	 */
	std::vector<PolyMesh::VertexHandle> face_vertices(const PolyMesh &mesh, const PolyMesh::FaceHandle fh) {
		vector<PolyMesh::VertexHandle> vhs;
		std::copy(mesh.cfv_ccwbegin(fh), mesh.cfv_ccwend(fh), vhs.begin());
		return vhs;
	}

	/**
	 * Get the points of a face.
	 * @param mesh The mesh.
	 * @param fh The face.
	 * @returns The points of the vertices of the given face.
	 */
	std::array<OpenMesh::Vec3d, 3> face_points(const TriMesh &mesh, const TriMesh::FaceHandle fh) {
		auto vhs = face_vertices(mesh, fh);
		array<OpenMesh::Vec3d, 3> points;
		std::transform(vhs.begin(), vhs.end(), points.begin(), [&](MishMesh::TriMesh::VertexHandle &vh) { return mesh.point(vh); });
		return points;
	}

	/**
	 * Get the points of a face.
	 * @param mesh The mesh.
	 * @param fh The face.
	 * @returns The points of the vertices of the given face.
	 */
	std::vector<OpenMesh::Vec3d> face_points(const PolyMesh &mesh, const PolyMesh::FaceHandle fh) {
		auto vhs = face_vertices(mesh, fh);
		vector<OpenMesh::Vec3d> points;
		std::transform(vhs.begin(), vhs.end(), std::back_inserter(points), [&](MishMesh::PolyMesh::VertexHandle &vh) { return mesh.point(vh); });
		return points;
	}

	/**
	 * Find the halfedge from vh1 to vh2.
	 *
	 * @param mesh The mesh in which to search for the halfedge.
	 * @param vh1 The first vertex handle.
	 * @param vh2 The second vertex handle.
	 * @returns The halfedge from vh1 to vh2 if found, otherwise an invalid handle.
	 */
	template<typename MeshT>
	OpenMesh::SmartHalfedgeHandle find_halfedge(const MeshT &mesh, const typename MeshT::VertexHandle vh1,
	                                            const typename MeshT::VertexHandle &vh2) {
		if(!vh1.is_valid()) {
			return OpenMesh::SmartHalfedgeHandle();
		}
		FOR_CVOH(h_it, vh1) {
			if(h_it->to() == vh2) {
				return *h_it;
			}
		}
		// If not found, return an invalid halfedge
		return OpenMesh::SmartHalfedgeHandle();
	}

	/**
	 * Find the edge connecting vh1 and vh2.
	 *
	 * @param mesh The mesh in which to search for the halfedge.
	 * @param vh1 The first vertex handle.
	 * @param vh2 The second vertex handle.
	 * @returns The edge between vh1 and vh2 if found, otherwise an invalid handle.
	 */
	template<typename MeshT>
	OpenMesh::SmartEdgeHandle find_edge(const MeshT &mesh, const typename MeshT::VertexHandle vh1,
	                                    const typename MeshT::VertexHandle &vh2) {
		const auto heh = find_halfedge(mesh, vh1, vh2);
		if(heh.is_valid()) {
			return heh.edge();
		}
		// If not found, return an invalid edge
		return OpenMesh::SmartEdgeHandle();
	}

	/**
	 * Split the edge belonging to a given halfedge handle and return the two halfedges
	 * ordered such that the to-vertex of the first edge is the from-vertex of the second.
	 *
	 * @param mesh The mesh.
	 * @param heh The halfedge handle of the edge to split.
	 * @param vh The vertex handle to insert into the edge.
	 * @return A pair of halfedge handles representing the two segments of the split edge.
	 * @note The input halfedge handle stays valid and is one of the two returned halfedge handles.
	 */
	SmartHalfedgePair split_edge(MishMesh::TriMesh &mesh, const OpenMesh::SmartHalfedgeHandle heh,
	                                                     const MishMesh::TriMesh::VertexHandle &vh) {
		auto heh0 = heh.edge().h0();
		mesh.split_edge(heh.edge(), vh);
		// The edge heh0 is split into: new_edge.h1 -> vh -> heh0
		// TriMesh::split_edge also adds two edges that re-triangulate the adjacent faces.
		// We find the new edge by circulating around the vertex and skipping the re-triangulation edge.
		auto new_heh = heh0.prev().opp().prev();
		if(heh == heh0) {
			return {new_heh, heh};
		} else {
			return {heh, new_heh.opp()};
		}
	}

	/**
	 * Split the edge belonging to a given halfedge handle and return the two halfedges
	 * ordered such that the to-vertex of the first edge is the from-vertex of the second.
	 *
	 * @param mesh The mesh.
	 * @param heh The halfedge handle of the edge to split.
	 * @param vh The vertex handle to insert into the edge.
	 * @return A pair of halfedge handles representing the two segments of the split edge.
	 * @note The input halfedge handle stays valid and is one of the two returned halfedge handles.
	 */
	MishMesh::SmartHalfedgePair split_edge(MishMesh::PolyMesh &mesh, const OpenMesh::SmartHalfedgeHandle heh,
	                                                      const MishMesh::PolyMesh::VertexHandle &vh) {
		auto eh = mesh.edge_handle(heh);
		auto heh0 = heh.edge().h0();
		mesh.split_edge(heh.edge(), vh);
		// heh0 is split into: new_edge.h0 -> vh -> heh0
		auto new_heh = heh0.prev();
		if(heh == heh0) {
			return {new_heh, heh};
		} else {
			return {heh, new_heh.opp()};
		}
	}

	/**
	 * Get the nearest vertex handle of a face from given barycentric coordinates
	 * @param mesh The mesh.
	 * @param fh The face handle.
	 * @param barycentric_coordinates Barycentric coordinates that sum to 1.0.
	 * @returns The nearest vertex to the point inside the triangle.
	 */
	template<typename VectorT>
	MishMesh::TriMesh::VertexHandle nearest_vh(const TriMesh &mesh, const MishMesh::TriMesh::FaceHandle fh, const VectorT barycentric_coordinates) {
		assert(barycentric_coordinates[0] + barycentric_coordinates[1] + barycentric_coordinates[2] == 1.0);

		auto vhs = face_vertices(mesh, fh);
		if(barycentric_coordinates[0] >= 1 / 3.) {
			return vhs[0];
		} else if(barycentric_coordinates[1] > 1 / 3.) {
			return vhs[1];
		} else {
			return vhs[2];
		}
	}

	/**
	 * Compute the triangle area from three points.
	 * @param points The coordinates of the triangle vertices.
	 * @tparam DIM The dimension of the ambiant space for 2D / 3D vectors.
	 * @returns the area of the triangle.
	 */
	template<int DIM>
	double compute_area(const array<OpenMesh::VectorT<double, DIM>, 3> points) {
		double l01 = (points[1] - points[0]).norm();
		double l12 = (points[2] - points[1]).norm();
		double l20 = (points[0] - points[2]).norm();
		double s = (l01 + l12 + l20) / 2.0;
		return sqrt(s * (s - l01) * (s - l12) * (s - l20));
	}

	/**
	 * Compute the triangle area from three vertex handles.
	 * @param mesh The mesh.
	 * @param vertices The triangle vertices.
	 * @returns the area of the triangle.
	 */
	double compute_area(const TriMesh &mesh, const array<TriMesh::VertexHandle, 3> vertices) {
		array<OpenMesh::Vec3d, 3> points;
		std::transform(vertices.begin(), vertices.end(), points.begin(), [&](const MishMesh::TriMesh::VertexHandle &vh) { return mesh.point(vh); });
		return compute_area(points);
	}

	/**
	 * Compute the triangle area of a face.
	 * @param mesh The mesh.
	 * @param fh The face handle.
	 * @returns the area of the triangle.
	 */
	double compute_area(const TriMesh &mesh, const TriMesh::FaceHandle fh) {
		return compute_area(mesh, face_vertices(mesh, fh));
	}

	/**
	 * Get the vertex with an obtuse angle in a triangle, if the face has an obtuse angle.
	 * @param mesh The mesh.
	 * @param fh The face handle.
	 * @returns A vertex handle to the vertex that belongs to the obtuse angle,
	 *          or an invalid vertex handle if the triangle has no obtuse angles.
	 */
	TriMesh::VertexHandle obtuse_vertex(const TriMesh &mesh, const TriMesh::FaceHandle fh) {
		std::array<double, 3> sqr_lengths;
		std::array<TriMesh::VertexHandle, 3> vhs;
		int i = 0;
		FOR_CFH(h_it, fh) {
			sqr_lengths[i] = mesh.calc_edge_sqr_length(*h_it);
			vhs[i] = mesh.opposite_vh(*h_it);
			i++;
		}
		size_t max_idx = std::distance(sqr_lengths.begin(), std::max_element(sqr_lengths.begin(), sqr_lengths.end()));
		if(sqr_lengths[max_idx] / (sqr_lengths[(max_idx + 1) % 3] + sqr_lengths[(max_idx + 2) % 3]) > 1.0) {
			return vhs[max_idx];
		}
		return TriMesh::InvalidVertexHandle;
	}

	template<typename MeshT>
	Flatness is_flat(const MeshT &mesh, double eps) {
		if(mesh.n_vertices() < 3) {
			return Flatness::DEGENERATE;
		}
		const auto p1 = mesh.point(mesh.vertex_handle(0));

		bool flat_x = true;
		bool flat_y = true;
		bool flat_z = true;
		for(const auto vh : mesh.vertices()) {
			const auto p2 = mesh.point(vh);
			if(abs(p1[0] - p2[0]) > eps) flat_x = false;
			if(abs(p1[1] - p2[1]) > eps) flat_y = false;
			if(abs(p1[2] - p2[2]) > eps) flat_z = false;

			if(flat_x == false && flat_y == false && flat_z == false) return Flatness::NONFLAT;
		}
		// When two dimensions are constant for all points, the points are collinear
		if(flat_x == true && flat_y == true) return Flatness::DEGENERATE;
		if(flat_y == true && flat_z == true) return Flatness::DEGENERATE;
		if(flat_z == true && flat_x == true) return Flatness::DEGENERATE;

		if(flat_x) {
			assert(!flat_y && !flat_z);
			return Flatness::X;
		}
		if(flat_y) {
			assert(!flat_x && !flat_z);
			return Flatness::Y;
		}
		if(flat_z) {
			assert(!flat_x && !flat_y);
			return Flatness::Z;
		}
		return Flatness::DEGENERATE;
	}

#ifdef HAS_EIGEN
	/**
	 * Convert a TriMesh to face and vertex vectors.
	 * @param[in] mesh The mesh.
	 * @returns A pair (V, F) with a Vx3 matrix of the vertices and a Fx3 matrix mapping faces to their vertices.
	 * @note The function assumes that the mesh indices are consistent. Call mesh.garbage_collection() before
	 *       using this function, when you elements from the mesh.
	 */
	template<typename VertexMatrixT, typename FaceMatrixT>
	std::pair<VertexMatrixT, FaceMatrixT> convert_to_face_vertex_mesh(const MishMesh::TriMesh &mesh) {
		VertexMatrixT V(mesh.n_vertices(), 3);
		for(int i = 0; i < V.rows(); i++) {
			V.row(i) = Eigen::Vector3d(mesh.point(mesh.vertex_handle(i)).data());
		}

		FaceMatrixT F(mesh.n_faces(), 3);
		for(int i = 0; i < F.rows(); i++) {
			const auto fh = mesh.face_handle(i);
			int v_idx = 0;
			FOR_CFV(v_it, fh) {
				F(i, v_idx) = v_it->idx();
				v_idx++;
			}
		}

		return std::make_pair(V, F);
	}

	template std::pair<Eigen::MatrixX3d, Eigen::MatrixX3i> convert_to_face_vertex_mesh(const MishMesh::TriMesh &mesh);
	template std::pair<Eigen::MatrixXd, Eigen::MatrixXi> convert_to_face_vertex_mesh(const MishMesh::TriMesh &mesh);
	template TriMesh::VertexHandle nearest_vh(const TriMesh &mesh, const MishMesh::TriMesh::FaceHandle fh, const Eigen::Vector3d barycentric_coordinates);
	template TriMesh::VertexHandle nearest_vh(const TriMesh &mesh, const MishMesh::TriMesh::FaceHandle fh, const Eigen::Vector3f barycentric_coordinates);
#endif

	template double compute_area(const std::array<OpenMesh::VectorT<double, 2>, 3> points);
	template double compute_area(const std::array<OpenMesh::VectorT<double, 3>, 3> points);
	template TriMesh::VertexHandle nearest_vh(const TriMesh &mesh, const MishMesh::TriMesh::FaceHandle fh, const std::array<double, 3> barycentric_coordinates);

	template Flatness is_flat(const TriMesh &mesh, double eps);
	template Flatness is_flat(const PolyMesh &mesh, double eps);

	// clang-format off
	template std::array<typename MishMesh::TriMesh::VertexHandle, 2> edge_vertices(const MishMesh::TriMesh &mesh, const typename MishMesh::TriMesh::HalfedgeHandle heh);
	template std::array<typename MishMesh::TriMesh::VertexHandle, 2> edge_vertices(const MishMesh::TriMesh &mesh, const typename MishMesh::TriMesh::EdgeHandle eh);
	template std::array<typename MishMesh::TriMesh::Point, 2> edge_points(const MishMesh::TriMesh &mesh, const typename MishMesh::TriMesh::HalfedgeHandle heh);
	template std::array<typename MishMesh::TriMesh::Point, 2> edge_points(const MishMesh::TriMesh &mesh, const typename MishMesh::TriMesh::EdgeHandle eh);

	template std::array<typename MishMesh::PolyMesh::VertexHandle, 2> edge_vertices(const MishMesh::PolyMesh &mesh, const typename MishMesh::PolyMesh::HalfedgeHandle heh);
	template std::array<typename MishMesh::PolyMesh::VertexHandle, 2> edge_vertices(const MishMesh::PolyMesh &mesh, const typename MishMesh::PolyMesh::EdgeHandle eh);
	template std::array<typename MishMesh::PolyMesh::Point, 2> edge_points(const MishMesh::PolyMesh &mesh, const typename MishMesh::PolyMesh::HalfedgeHandle heh);
	template std::array<typename MishMesh::PolyMesh::Point, 2> edge_points(const MishMesh::PolyMesh &mesh, const typename MishMesh::PolyMesh::EdgeHandle eh);

	template OpenMesh::SmartHalfedgeHandle find_halfedge(const MishMesh::TriMesh &mesh, const typename MishMesh::TriMesh::VertexHandle vh1, const typename MishMesh::TriMesh::VertexHandle &vh2);
	template OpenMesh::SmartHalfedgeHandle find_halfedge(const MishMesh::PolyMesh &mesh, const typename MishMesh::PolyMesh::VertexHandle vh1, const typename MishMesh::PolyMesh::VertexHandle &vh2);
	template OpenMesh::SmartEdgeHandle find_edge(const MishMesh::TriMesh &mesh, const typename MishMesh::TriMesh::VertexHandle vh1, const typename MishMesh::TriMesh::VertexHandle &vh2);
	template OpenMesh::SmartEdgeHandle find_edge(const MishMesh::PolyMesh &mesh, const typename MishMesh::PolyMesh::VertexHandle vh1, const typename MishMesh::PolyMesh::VertexHandle &vh2);
}
