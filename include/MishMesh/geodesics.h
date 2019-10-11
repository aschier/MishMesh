#pragma once

#include <MishMesh/TriMesh.h>
#include <exception>

namespace MishMesh {
	class NoOverlap: public std::exception {};
	typedef OpenMesh::VPropHandleT<double> GeodesicDistanceProperty;
	enum GeodesicType {
		NOVOTNI,
		HEAT,
	};

	std::pair<OpenMesh::Vec2d, OpenMesh::Vec2d> compute_projected_origins(double edge_length, double T1, double T2);
	template<int DIM>
	double compute_distance(const OpenMesh::VectorT<double, DIM> p, const OpenMesh::VectorT<double, DIM> p1, const OpenMesh::VectorT<double, DIM> p2, const double T1, const double T2);
	double compute_distance(const TriMesh &mesh, const TriMesh::VertexHandle &vh3, const TriMesh::VertexHandle &edge_vh1, const TriMesh::VertexHandle &edge_vh2, const OpenMesh::VPropHandleT<double> &distProp);

	void compute_novotni_geodesics(TriMesh &mesh, const TriMesh::VertexHandle start_vh, const GeodesicDistanceProperty geodesicDistanceProperty);
#ifdef HAS_EIGEN
	void compute_heat_geodesics(TriMesh &mesh, const TriMesh::VertexHandle start_vh, GeodesicDistanceProperty geodesicDistanceProperty, double t = 0.1);
#endif

	inline void compute_geodesics(TriMesh &mesh, const TriMesh::VertexHandle start_vh, GeodesicDistanceProperty geodesicDistanceProperty, GeodesicType geodesicType = NOVOTNI) {
		if(geodesicType == NOVOTNI) {
			compute_novotni_geodesics(mesh, start_vh, geodesicDistanceProperty);
#ifdef HAS_EIGEN
		} else if(geodesicType == HEAT) {
			compute_heat_geodesics(mesh, start_vh, geodesicDistanceProperty);
#endif
		} else {
			assert(false);
		}
	}
}
