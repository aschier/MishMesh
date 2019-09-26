#pragma once

#include <queue>
#include <set>

#include <MishMesh/PolyMesh.h>
#include <MishMesh/TriMesh.h>
#include <MishMesh/OpenMesh/DoublePrecisionTraits.h>

namespace MishMesh {
	/**
	 * Stores the result of a dijkstra shortest path search.
	 */
	template<typename MeshT>
	struct DijkstraResult {
		/// length The lenght of the path or std::numeric_limits<double>::infinity() when no path was found.
		double length;
		/// An ordered list of edge handles describing the path.
		std::vector<typename MeshT::EdgeHandle> edges;
		/// An ordered list of vertex handles describing the path.
		std::vector<typename MeshT::VertexHandle> vertices;
		/// Returns true, if a path was found.
		bool is_valid() {
			return length < std::numeric_limits<double>::infinity() && length > 0;
		}
	};

	/**
	 * The edges between two vertices.
	 */
	template<typename MeshT>
	struct EdgePath {
		int v1, v2;
		double length;
		std::vector<typename MeshT::EdgeHandle> edges;
	};

	/**
	 * An edge on a path. PathEdges can be ordered by the length of the path up to this edge.
	 * @param mesh The mesh the edge belongs to. Used for length calculations.
	 * @param The property handle for the "shortest path" property on the mesh.
	 * @param halfedge_handle The handle to the edge in the mesh.
	 */
	template<typename MeshT>
	class PathEdge {
	public:
		MeshT *mesh;
		OpenMesh::HPropHandleT<double> *prop_edge_shortest_path_length;
		typename MeshT::HalfedgeHandle halfedge_handle;
	};
	/**
	 * A comparator, that compares two PathEdges by the length of the path up to this edge.
	 */
	template<typename MeshT>
	class GreaterPathlengh {
	public:
		bool operator()(const PathEdge<MeshT> &a, const PathEdge<MeshT> &b) const {
			return a.mesh->property(*a.prop_edge_shortest_path_length, a.halfedge_handle) > b.mesh->property(*b.prop_edge_shortest_path_length, b.halfedge_handle);
		}
	};


	/**
	 * Calculate the euclidean distance between the two endpoints of the edge.
	 */
	template<typename MeshT>
	inline double edge_length(const MeshT &mesh, const typename MeshT::HalfedgeHandle edge, const void *param) {
		return mesh.calc_edge_length(edge);
	}

	template<typename MeshT>
	DijkstraResult<MeshT> dijkstra(const typename MeshT::VertexHandle start_vh, const typename MeshT::VertexHandle target_vh, MeshT &mesh, double edge_cost_function(MeshT &mesh, const typename MeshT::HalfedgeHandle edge, const void *param) = edge_length<MeshT>, void *edge_cost_param = nullptr);
}
