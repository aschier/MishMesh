#pragma once

#include <MishMesh/BBox.h>
#include <MishMesh/PolyMesh.h>
#include <functional>

namespace MishMesh {
	namespace Meshing {
		/// A signed distance function, that returns the distance and the normal vector for a given point.
		typedef std::function<std::pair<double, OpenMesh::Vec3d>(const OpenMesh::Vec3d, void *)> DistanceFunction;

		/**
		 * The method how vertices should be fit inside a voxel.
		 * - None: Use the center of the voxel. This is usually not useful.
		 * - Average: Average the distances along the edges to determine the position. This option gives smoother meshes.
		 * - SVD: Fit the vertex to a plane defined by the normals at the intersection points. This option allows to preserve sharp edges.
		 */
#ifdef HAS_EIGEN
		enum VertexFit {
			None,
			Average,
			SVD
		};
#else
		enum VertexFit {
			None,
			Average
		};
#endif

		/**
		* Create a mesh using a signed distance function that returns hermite data (the distance and a normal) for a point.
		 * @param distance_function The distance function.
		 * @param distance_function_data optional data for the distance function or a nullptr.
		 * @param mesh_bounding_box The bounding box of the mesh.
		 * @param resolution The resolution for a uniform grid inside the bounding box.
		 * @param vertexFit The method to fit the vertices.
		 * @returns a quad mesh of the signed distance function.
		 */
		MishMesh::PolyMesh dual_contouring(DistanceFunction distance_function, void *distance_function_data = nullptr, MishMesh::BBox<OpenMesh::Vec3d, 3> mesh_bounding_box = {{-1.0, -1.0, -1.0}, {1.0, 1.0, 1.0}}, OpenMesh::Vec3i resolution = {50, 50, 50}, VertexFit vertexFit = VertexFit::Average);
	}
}