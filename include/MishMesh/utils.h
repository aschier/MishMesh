#pragma once
#include <MishMesh/TriMesh.h>
#include <MishMesh/PolyMesh.h>
#include <array>
#include <limits>

#ifdef HAS_EIGEN
#include <Eigen/Eigen>
#endif

namespace MishMesh {
	TriMesh::HalfedgeHandle opposite_halfedge(const TriMesh &mesh, const TriMesh::FaceHandle &fh, const TriMesh::VertexHandle &vh);

	std::array<TriMesh::VertexHandle, 2> edge_vertices(const TriMesh &mesh, const TriMesh::HalfedgeHandle heh);
	std::array<TriMesh::VertexHandle, 2> edge_vertices(const TriMesh &mesh, const TriMesh::EdgeHandle eh);

	std::array<OpenMesh::Vec3d, 2> edge_points(const TriMesh &mesh, const TriMesh::HalfedgeHandle heh);
	std::array<OpenMesh::Vec3d, 2>  edge_points(const TriMesh &mesh, const TriMesh::EdgeHandle eh);

	std::array<TriMesh::VertexHandle, 3> face_vertices(const TriMesh &mesh, const TriMesh::FaceHandle fh);
	std::vector<PolyMesh::VertexHandle> face_vertices(const PolyMesh &mesh, const PolyMesh::FaceHandle fh);

	std::array<OpenMesh::Vec3d, 3> face_points(const TriMesh &mesh, const TriMesh::FaceHandle fh);
	std::vector<OpenMesh::Vec3d> face_points(const PolyMesh &mesh, const PolyMesh::FaceHandle fh);

	template<typename VectorT>
	TriMesh::VertexHandle nearest_vh(const TriMesh &mesh, const MishMesh::TriMesh::FaceHandle fh, const VectorT barycentric_coordinates);

	template<int DIM>
	double compute_area(const std::array<OpenMesh::VectorT<double, DIM>, 3> points);
	double compute_area(const TriMesh &mesh, const std::array<TriMesh::VertexHandle, 3> vertices);
	double compute_area(const TriMesh &mesh, const TriMesh::FaceHandle fh);

	TriMesh::VertexHandle obtuse_vertex(const TriMesh &mesh, const TriMesh::FaceHandle fh);

	enum Flatness {
		X,
		Y,
		Z,
		NONFLAT,
		DEGENERATE
	};

	/** Test if a mesh is flat.
	 * @param mesh The mesh.
	 * @param eps An optional epsilon for ignoring small differences
	 * @tparam MeshT the mesh type
	 * @returns Flatness::X, Flatness:Y, Flatness::Z if the mesh is flat in one direction,
	 *          Flatness::NONFLAT if the mesh is not flat and Flatness::DEGENERATE when the mesh
	 *          is flat in more than one direction, i.e., all points are collinar, or the mesh has
	 *          no faces.
	 */
	template<typename MeshT>
	Flatness is_flat(const MeshT &mesh, double eps = 0);

	/**
	 * Test if a face has an obtuse angle.
	 * @param mesh The mesh.
	 * @param fh a face handle in the mesh.
	 * @returns True, if the face has an obtuse angle.
	 */
	inline bool is_obtuse(const TriMesh &mesh, const TriMesh::FaceHandle fh) {
		return obtuse_vertex(mesh, fh).is_valid();
	};

	/**
	 * Test if a vertex of a face belongs to an obtuse angle.
	 * @param mesh The mesh.
	 * @param fh a face handle in the mesh.
	 * @param vh a vertex handle in the mesh.
	 * @returns True, if the vertex belongs to an obtuse angle in the face.
	 * @note The method does not test if the vertex belongs to the face and
	 *       always returns false when the vertex does not belong to the face.
	 */
	inline bool is_obtuse(const TriMesh &mesh, const TriMesh::FaceHandle fh, const TriMesh::VertexHandle vh) {
		return obtuse_vertex(mesh, fh) == vh;
	};

	/**
	 * Computes the euler characteristic of a mesh.
	 * @param mesh The mesh.
	 * @returns The euler characteristic.
	 */
	template<typename MeshT>
	inline size_t euler_characteristic(const MeshT &mesh) {
		return mesh.n_vertices() - mesh.n_edges() + mesh.n_faces();
	}

#ifdef HAS_EIGEN
	template<typename VertexMatrixT = Eigen::MatrixX3d, typename FaceMatrixT = Eigen::MatrixX3i>
	std::pair<VertexMatrixT, FaceMatrixT> convert_to_face_vertex_mesh(const MishMesh::TriMesh &mesh);
#endif
}
