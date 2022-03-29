#pragma once

#include <queue>
#include <set>

#include <MishMesh/PolyMesh.h>
#include <MishMesh/TriMesh.h>

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
		PathEdge(MeshT *mesh,
		         OpenMesh::HPropHandleT<double> *prop_edge_shortest_path_length,
		         typename MeshT::HalfedgeHandle halfedge_handle) : mesh(mesh), prop_edge_shortest_path_length(prop_edge_shortest_path_length), halfedge_handle(halfedge_handle),
		                                                           start_vh(&MeshT::InvalidVertexHandle), target_vh(&MeshT::InvalidVertexHandle){};
		PathEdge(MeshT *mesh,
		         OpenMesh::HPropHandleT<double> *prop_edge_shortest_path_length,
		         typename MeshT::HalfedgeHandle halfedge_handle,
		         const typename MeshT::VertexHandle *start_vh,
		         const typename MeshT::VertexHandle *target_vh) : mesh(mesh), prop_edge_shortest_path_length(prop_edge_shortest_path_length), halfedge_handle(halfedge_handle),
		                                                          start_vh(start_vh), target_vh(target_vh){};
		MeshT *mesh;
		const MishMesh::TriMesh::VertexHandle *start_vh;
		const MishMesh::TriMesh::VertexHandle *target_vh;
		OpenMesh::HPropHandleT<double> *prop_edge_shortest_path_length;
		typename MeshT::HalfedgeHandle halfedge_handle;
	};

	/**
	 * A comparator, that compares two PathEdges by the length of the path up to this edge.
	 */
	template<typename MeshT = MishMesh::TriMesh>
	class GreaterPathlength {
	public:
		bool operator()(const PathEdge<MeshT> &a, const PathEdge<MeshT> &b) const {
			return a.mesh->property(*a.prop_edge_shortest_path_length, a.halfedge_handle) > b.mesh->property(*b.prop_edge_shortest_path_length, b.halfedge_handle);
		}
	};

	/**
	 * A comparator, that compares two PathEdges by the length of the path up to this edge plus the l2 distance to a target vertex
	 */
	template<typename MeshT = MishMesh::TriMesh>
	class L2HeuristicComparator {
	public:
		bool operator()(const PathEdge<MeshT> &a, const PathEdge<MeshT> &b) const {
			assert(a.target_vh == b.target_vh);
			double lengthA = a.mesh->property(*a.prop_edge_shortest_path_length, a.halfedge_handle);
			lengthA += (a.mesh->point(*a.target_vh) - a.mesh->point(a.mesh->to_vertex_handle(a.halfedge_handle))).norm();
			double lengthB = b.mesh->property(*b.prop_edge_shortest_path_length, b.halfedge_handle);
			lengthB += (b.mesh->point(*b.target_vh) - b.mesh->point(b.mesh->to_vertex_handle(b.halfedge_handle))).norm();
			return lengthA > lengthB;
		}
	};

	/**
	 * Calculate the euclidean distance between the two endpoints of the edge.
	 */
	template<typename MeshT>
	inline double edge_length(const MeshT &mesh, const typename MeshT::HalfedgeHandle edge, const void *param) {
		return mesh.calc_edge_length(edge);
	}

	/**
	 * Calculate the euclidean distance between the two endpoints of the edge.
	 */
	template<typename MeshT>
	inline double edge_length(MeshT &mesh, const typename MeshT::HalfedgeHandle edge, const void *param) {
		return mesh.calc_edge_length(edge);
	}

	template<typename MeshT, typename ComparatorT = GreaterPathlength<MeshT>>
	DijkstraResult<MeshT> dijkstra(const typename MeshT::VertexHandle start_vh, const typename MeshT::VertexHandle target_vh, MeshT &mesh, double edge_cost_function(MeshT &mesh, const typename MeshT::HalfedgeHandle edge, const void *param) = edge_length<MeshT>, void *edge_cost_param = nullptr);

	template<typename MeshT>
	DijkstraResult<MeshT> trace_path(const MeshT &mesh,
	                          const typename MeshT::VertexHandle &target_vh,
	                          const OpenMesh::VPropHandleT<double> &prop_vertex_shortest_path_length,
	                          const OpenMesh::HPropHandleT<double> &prop_edge_shortest_path_length);

	};

	/**
	 * Calculate the euclidean distance between the two endpoints of the edge.
	 */
	template<typename MeshT>
	inline double edge_length(const MeshT &mesh, const typename MeshT::HalfedgeHandle edge, const void *param) {
		return mesh.calc_edge_length(edge);
	}

	/**
	 * Calculate the euclidean distance between the two endpoints of the edge.
	 */
	template<typename MeshT>
	inline double edge_length(MeshT &mesh, const typename MeshT::HalfedgeHandle edge, const void *param) {
		return mesh.calc_edge_length(edge);
	}
