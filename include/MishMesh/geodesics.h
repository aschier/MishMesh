#pragma once

#include <MishMesh/TriMesh.h>
#include <exception>

namespace MishMesh {
	class NoOverlap: public std::exception {};
	typedef OpenMesh::VPropHandleT<double> GeodesicDistanceProperty;

	void compute_geodesics(TriMesh &mesh, const TriMesh::VertexHandle start_vh, const GeodesicDistanceProperty geodesicDistanceProperty);
	std::pair<OpenMesh::Vec2d, OpenMesh::Vec2d> compute_projected_origins(double edge_length, double T1, double T2);

	template<int DIM>
	double compute_distance(const OpenMesh::VectorT<double, DIM> p, const OpenMesh::VectorT<double, DIM> p1, const OpenMesh::VectorT<double, DIM> p2, const double T1, const double T2);
	double compute_distance(const TriMesh &mesh, const TriMesh::VertexHandle &vh3, const TriMesh::VertexHandle &edge_vh1, const TriMesh::VertexHandle &edge_vh2, const OpenMesh::VPropHandleT<double> &distProp);
}

