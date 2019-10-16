#include "MishMesh/utils.h"
#include <MishMesh/macros.h>

#include <MishMesh/TriMesh.h>
#include <MishMesh/PolyMesh.h>

using namespace std;

namespace MishMesh {
	/**
	 * Get the halfedge of a triangle, that is opposite to a given vertex.
	 * @param mesh The mesh.
	 * @param fh The face.
	 * @param vh The vertex.
	 * @returns The halfedge of the face, that is opposite to the vertex.
	 */
	TriMesh::HalfedgeHandle opposite_halfedge(const TriMesh &mesh, const TriMesh::FaceHandle &fh, TriMesh::VertexHandle &vh) {
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
			vhs[j++] = *v_it;
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
		std::transform(vhs.begin(), vhs.end(), points.begin(), [&](MishMesh::TriMesh::VertexHandle &vh) {return mesh.point(vh); });
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
		std::transform(vhs.begin(), vhs.end(), points.begin(), [&](MishMesh::TriMesh::VertexHandle &vh) {return mesh.point(vh); });
		return points;
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
		return sqrt(s*(s - l01)*(s - l12)*(s - l20));
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
		for(auto v : vertices){
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

	template double compute_area(const std::array<OpenMesh::VectorT<double, 2>, 3> points);
	template double compute_area(const std::array<OpenMesh::VectorT<double, 3>, 3> points);
}
