#pragma once

#include <MishMesh/TriMesh.h>
#include <exception>
#include <OpenMesh/Tools/Utils/HeapT.hh>

namespace MishMesh {
	class NoOverlap: public std::exception {};
	typedef OpenMesh::VPropHandleT<double> GeodesicDistanceProperty;
	typedef OpenMesh::VPropHandleT<int> HeapIndexProperty;

	template<typename HeapEntry>
	struct VertexHeapInterface {
		 VertexHeapInterface(MishMesh::TriMesh &_mesh, const GeodesicDistanceProperty _propDistance, const HeapIndexProperty _propHeapIndex):
			  mesh(_mesh), propDistance(_propDistance), propHeapIndex(_propHeapIndex) {
			  for(auto &vh : mesh.vertices()) {
				   mesh.property(propHeapIndex, vh) = -1;
			  }
		 }
		 /// Comparison of two HeapEntry's: strict less
		 bool less(const HeapEntry& _e1, const HeapEntry& _e2) {
			  return mesh.property(propDistance, _e1) < mesh.property(propDistance, _e2);
		 };

		 /// Comparison of two HeapEntry's: strict greater
		 bool greater(const HeapEntry& _e1, const HeapEntry& _e2) {
			  return mesh.property(propDistance, _e1) > mesh.property(propDistance, _e2);
		 }

		 /// Get the heap position of HeapEntry _e
		 int get_heap_position(const HeapEntry& _e) {
			  return mesh.property(propHeapIndex, _e);
		 };

		 /// Set the heap position of HeapEntry _e
		 void set_heap_position(HeapEntry& _e, int _i) {
			  mesh.property(propHeapIndex, _e) = _i;
		 };

		 const GeodesicDistanceProperty propDistance;
		 const HeapIndexProperty propHeapIndex;
		 TriMesh &mesh;
	};
	typedef OpenMesh::Utils::HeapT<MishMesh::TriMesh::VertexHandle, VertexHeapInterface<MishMesh::TriMesh::VertexHandle>> VertexHeap;

	inline void insert_or_update(VertexHeap heap, TriMesh::VertexHandle vh) {
		 if(heap.is_stored(vh)) {
			  heap.update(vh);
		 } else{
			  heap.insert(vh);
		 };
	}


	enum GeodesicType {
		NOVOTNI,
		HEAT,
	};

	std::pair<OpenMesh::Vec2d, OpenMesh::Vec2d> compute_projected_origins(double edge_length, double T1, double T2);
	template<int DIM>
	double compute_distance(const OpenMesh::VectorT<double, DIM> p, const OpenMesh::VectorT<double, DIM> p1, const OpenMesh::VectorT<double, DIM> p2, const double T1, const double T2);
	double compute_distance(const TriMesh &mesh, const TriMesh::VertexHandle &vh3, const TriMesh::VertexHandle &edge_vh1, const TriMesh::VertexHandle &edge_vh2, const OpenMesh::VPropHandleT<double> &distProp);

	void compute_novotni_geodesics(TriMesh &mesh, const TriMesh::VertexHandle start_vh, const GeodesicDistanceProperty geodesicDistanceProperty, const bool handle_obtuse=true);
#ifdef HAS_EIGEN
	void compute_heat_geodesics(TriMesh &mesh, const TriMesh::VertexHandle start_vh, GeodesicDistanceProperty geodesicDistanceProperty, double t = 0.1);
	void compute_heat_geodesics(TriMesh &mesh, const std::vector<TriMesh::VertexHandle> start_vhs, GeodesicDistanceProperty geodesicDistanceProperty, double t = 0.1);
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
