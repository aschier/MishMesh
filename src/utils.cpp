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
	std::array<TriMesh::VertexHandle, 2> edge_vertices(const TriMesh &mesh, const TriMesh::HalfedgeHandle heh) {
		return {mesh.from_vertex_handle(heh), mesh.to_vertex_handle(heh)};
	}

	/**
	 * Get the vertex handles of an edge
	 * @param mesh The mesh.
	 * @param eh The edge.
	 * @returns The vertex handles of the vertices of the given edge.
	 */
	std::array<TriMesh::VertexHandle, 2> edge_vertices(const TriMesh &mesh, const TriMesh::EdgeHandle eh) {
		return edge_vertices(mesh, mesh.halfedge_handle(eh, 0));
	}

	/**
	 * Get the points of a half edge
	 * @param mesh The mesh.
	 * @param heh The halfedge.
	 * @returns The points of the vertices of the given halfedge.
	 */
	std::array<OpenMesh::Vec3d, 2> edge_points(const TriMesh &mesh, const TriMesh::HalfedgeHandle heh) {
		return {mesh.point(mesh.from_vertex_handle(heh)), mesh.point(mesh.to_vertex_handle(heh))};
	}

	/**
	 * Get the points of an edge
	 * @param mesh The mesh.
	 * @param eh The edge.
	 * @returns The points of the vertices of the given halfedge.
	 */
	std::array<OpenMesh::Vec3d, 2> edge_points(const TriMesh &mesh, const TriMesh::EdgeHandle eh) {
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
		short j = 0;
		FOR_CFV(v_it, fh) {
			vhs[j++] = *v_it;
		}
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
		short j = 0;
		FOR_CFV(v_it, fh) {
			vhs.push_back(*v_it);
		}
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
		std::transform(vhs.begin(), vhs.end(), points.begin(), [&](MishMesh::TriMesh::VertexHandle &vh) { return mesh.point(vh); });
		return points;
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
		short j = 0;
		for(auto v : vertices) {
			points[j++] = mesh.point(v);
		}
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
}
