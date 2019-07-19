#pragma once

#include <queue>
#include <set>

#include <OpenMesh/Core/Mesh/TriMeshT.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <MishMesh/OpenMesh/DoublePrecisionTraits.h>

namespace MishMesh {
	/// A triangle mesh structure with double precision values.
	typedef OpenMesh::TriMesh_ArrayKernelT<OpenMesh::DoublePrecisionTraits> TriMesh;

	/**
	 * Stores the result of a dijkstra shortest path search.
	 */
	struct DijkstraResult {
		/// length The lenght of the path or std::numeric_limits<double>::infinity() when no path was found.
		double length;
		/// An ordered list of edge handles describing the path.
		std::vector<TriMesh::EdgeHandle> edges;
		/// An ordered list of vertex handles describing the path.
		std::vector<TriMesh::VertexHandle> vertices;
		/// Returns true, if a path was found.
		bool is_valid() {
			return length < std::numeric_limits<double>::infinity() && length > 0;
		}
	};

	/**
	 * The edges between two vertices.
	 */
	struct EdgePath {
		int v1, v2;
		double length;
		std::vector<TriMesh::EdgeHandle> edges;
	};

	/**
	 * An edge on a path. PathEdges can be ordered by the length of the path up to this edge.
	 * @param mesh The mesh the edge belongs to. Used for length calculations.
	 * @param The property handle for the "shortest path" property on the mesh.
	 * @param halfedge_handle The handle to the edge in the mesh.
	 */
	class PathEdge {
	public:
		TriMesh *mesh;
		OpenMesh::HPropHandleT<double> *prop_edge_shortest_path_length;
		TriMesh::HalfedgeHandle halfedge_handle;
	};
	/**
	 * A comparator, that compares two PathEdges by the length of the path up to this edge.
	 */
	class GreaterPathlengh {
	public:
		bool operator()(const PathEdge &a, const PathEdge &b) const {
			return a.mesh->property(*a.prop_edge_shortest_path_length, a.halfedge_handle) > b.mesh->property(*b.prop_edge_shortest_path_length, b.halfedge_handle);
		}
	};


	/**
	 * Calculate the euclidean distance between the two endpoints of the edge.
	 */
	inline double edge_length(TriMesh &mesh, const TriMesh::HalfedgeHandle edge, const void *param) {
		return mesh.calc_edge_length(edge);
	}

	DijkstraResult dijkstra(const TriMesh::VertexHandle start_vh, const TriMesh::VertexHandle target_vh, TriMesh &mesh, double edge_cost_function(TriMesh &mesh, const TriMesh::HalfedgeHandle edge, const void *param) = edge_length, void *edge_cost_param = nullptr);
}
