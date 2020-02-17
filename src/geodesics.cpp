#include "MishMesh/geodesics.h"
#include <MishMesh/macros.h>
#include <MishMesh/transformations.h>
#include <MishMesh/utils.h>
#include <OpenMesh/Tools/Utils/HeapT.hh>

using namespace std;
using namespace MishMesh;


/**
 * Compute the two possible projected origin points for known distances from two points lying on the X-axis.
 * For symmetry reasons, there are two possible origins, that are (ox, oy) and (ox, -oy).
 *
 * @param edge_length The distance between the two points.
 * @param T1 The distance of the first point from the origin.
 * @param T2 The distance of the second point from the origin.
 * @returns a pair with the two possible origin points.
 * @throws NoOverlap when there are no points with the given distances from the two vertices, i.e,
 *         the two circles with radius T1, T2 around the vertices have no intersection.
 */
pair<OpenMesh::Vec2d, OpenMesh::Vec2d> MishMesh::compute_projected_origins(double edge_length_2, double T1_2, double T2_2) {
	OpenMesh::Vec2d o1;
	double A = 2 * T1_2 * edge_length_2 - pow(edge_length_2, 2) + 2 * T2_2 * edge_length_2;
	double B = pow(T1_2 - T2_2, 2);
	assert(isfinite(A));
	assert(isfinite(B));

	if(B > A) {
		throw NoOverlap();
	}

	double edge_length = sqrt(edge_length_2);
	double ox = 0.5 * (edge_length_2 + T1_2 - T2_2) / edge_length;
	double oy = 0.5 * sqrt(A - B) / edge_length;
	return make_pair(OpenMesh::Vec2d{ox, oy}, OpenMesh::Vec2d{ox, -oy});
}

/**
 * Compute the geodesic distance of the point p to the origin from points p1, p2 and their distances T1, T2 to the origin.
 * The method projects the origin into the plane of the triangle and then calculcates the new distance as the distance between
 * p and the projected origin.
 *
 * @param p The point for which the distance should be calculcated.
 * @param p1 A point with known distance.
 * @param p2 A second point with known distance.
 * @param T1 The distance between p1 and the origin.
 * @param T2 The distance between p2 and the origin.
 * @returns The geodesic distance of p to the origin.
 */
template<int DIM>
double MishMesh::compute_distance(const OpenMesh::VectorT<double, DIM> p, const OpenMesh::VectorT<double, DIM> p1, const OpenMesh::VectorT<double, DIM> p2, const double T1_2, const double T2_2) {
	auto points = embed_triangle(p1, p2, p);
	double v2x_2 = (p1 - p2).sqrnorm();
	auto v3 = points[2];

	assert(isfinite(T1_2));
	assert(isfinite(T2_2));

	try{
		OpenMesh::Vec2d o1, o2;
		std::tie(o1, o2) = compute_projected_origins(v2x_2, T1_2, T2_2);
		double T3_2 = max((v3 - o1).sqrnorm(), (v3 - o2).sqrnorm());
		assert(isfinite(T3_2));
		return T3_2;
	} catch(NoOverlap) {
		throw NoOverlap();
	}
}

/**
 * Compute the geodesic distance from the origin for a triangle vertex given by vh3, using the known geodesic distances of the vertices edge_vh1 and edge_vh2 of
 * the edge opposite to vh3.
 *
 * @param mesh The mesh.
 * @param vh3 The vertex for which the geodesic distance should be computed.
 * @param edge_vh1 The first vertex of the edge opposite to vh3.
 * @param edge_vh2 The second vertex of the edge opposite to vh3.
 * @param distProp A mesh property storing the geodesic distances from an origin.
 * @returns The geodesic distance of the vertex from the origin.
 * @note The distances from edge_vh1 and edge_vh2 must be already computed for this function to work.
 */
double MishMesh::compute_distance(const TriMesh &mesh, const TriMesh::VertexHandle &vh3, const TriMesh::VertexHandle &edge_vh1, const TriMesh::VertexHandle &edge_vh2, const OpenMesh::VPropHandleT<double> &distProp) {
	auto p1 = mesh.point(edge_vh1);
	auto p2 = mesh.point(edge_vh2);
	auto p = mesh.point(vh3);

	double T1_2 = mesh.property(distProp, edge_vh1);
	double T2_2 = mesh.property(distProp, edge_vh2);
	assert(isfinite(T1_2) && isfinite(T2_2));
	return compute_distance(p, p1, p2, T1_2, T2_2);
}

/**
 * Update the geodesic distance value of a vertex and reorder the close_vertices heap.
 * @param mesh The mesh.
 * @param update_vh The vertex handle to update.
 * @param new_distance The new distance value.
 * @param close_vertices A vertex heap, that is updated with the new distance.
 * @param geodesicDistanceProperty A vertex property storing the geodesic distance.
 * @note When the vertex heap does not contain update_vh, it is inserted.
 */
void update_distance(MishMesh::TriMesh &mesh, const MishMesh::TriMesh::VertexHandle update_vh, const double new_distance, VertexHeap &close_vertices, const GeodesicDistanceProperty geodesicDistanceProperty) {
	assert(!mesh.status(update_vh).tagged()); // The vertex is not fixed

	const double old_distance = mesh.property(geodesicDistanceProperty, update_vh); // current best distance
	if(new_distance > old_distance) return;

	mesh.property(geodesicDistanceProperty, update_vh) = new_distance;

	// If the vertex was unprocessed, add it to the close set.
	if(mesh.status(update_vh).tagged2()) {
		assert(!close_vertices.is_stored(update_vh));
		mesh.status(update_vh).set_tagged2(false);
		close_vertices.insert(update_vh);
		return;
	} else {
		assert(!mesh.status(update_vh).tagged());
		assert(close_vertices.is_stored(update_vh));
		close_vertices.update(update_vh);
		return;
	}
}

/**
 * Compute geodesic distances on a mesh using the Method of Novotni and Klein
 * [Novotni, M., & Klein, R. (2002). Computing geodesic distances on triangular meshes. In In Proc. of WSCG’2002.]
 * http://cg.cs.uni-bonn.de/de/publikationen/paper-details/novotni-2002-computing/
 *
 * @param[inout] mesh The mesh.
 * @param[in] start_vh A valid vertex in the mesh, that will be used as start vertex.
 * @param[in] geodesicDistanceProperty A mesh property to store the geodesic distances. The method assumes, that the property is already added to the mesh.
 */
void MishMesh::compute_novotni_geodesics(TriMesh &mesh, const TriMesh::VertexHandle start_vh, const GeodesicDistanceProperty geodesicDistanceProperty) {
	/*
	  The algorithm has three sets of vertices, fixed, close and unprocessed, and uses a region growing strategy similar to Dijkstra's algorithm to find
	  shortest geodesic distances. When Dijkstra's algorithm encounters a vertex, it uses the edge length between the vertex with known distance and the new vertex
	  to calculate the path length. This algorithm uses the known geodesic distances of two points adjacent to the new point, to project the origin into the plane
	  of the three points and then measures the distance from the new point to the projected origin.

	  The algorithm maintains three sets of vertices:
	  - A fixed vertex has a fixed distance, that will not change anymore.
	  - A close vertex is adajacent to a fixed vertex and already has a value, that still may change.
	  - Unprocessed vertices are not seen yet and will be added as close vertices, when a neighbor vertex is processed.
	  In each iteration, a trial vertex from the close set is added to the fixed set and all its neighbors are assigned the minimal distance between the shortest
	  current shortest path to the trial vertex plus the edge length from the trial vertex to the neighbor vertex.
	 */

	mesh.request_vertex_status();

	HeapIndexProperty propHeapIndex;
	mesh.add_property(propHeapIndex);

	VertexHeapInterface<MishMesh::TriMesh::VertexHandle> heapInterface(mesh, geodesicDistanceProperty, propHeapIndex);
	VertexHeap close_vertices(heapInterface);

	for(auto vh : mesh.vertices()) {
		mesh.status(vh).set_tagged(false); // True, when the vertex is fixed
		mesh.status(vh).set_tagged2(true); // True, when the vertex is unprocessed
		mesh.property(geodesicDistanceProperty, vh) = numeric_limits<double>::infinity();
	}

	// Assign the start vertex distance 0 and mark it as fixed
	mesh.status(start_vh).set_tagged(true);
	mesh.property(geodesicDistanceProperty, start_vh) = 0;

	// Set the distance of the neighbor vertices from the edge length and change them from unprocesses to close
	FOR_CVOH(h_it, start_vh) {
		auto vh = mesh.to_vertex_handle(*h_it);
		double distance = mesh.calc_edge_sqr_length(*h_it);
		mesh.property(geodesicDistanceProperty, vh) = distance;
		assert(!close_vertices.is_stored(vh));
		close_vertices.insert(vh);
		mesh.status(vh).set_tagged2(false);
	}

	// In each iteration, we pick a trial vertex from the close set and update the distances of its neighbor vertices, that are close or unprocessed
	while(!close_vertices.empty()) {
		auto trial_vh = close_vertices.front();
		close_vertices.pop_front();
		assert(isfinite(mesh.property(geodesicDistanceProperty, trial_vh)));
		mesh.status(trial_vh).set_tagged(true);
		// For each face adjacent to the trial vertex that contains a close or unprocessed vertex and a fixed vertex,
		// update the distance value of the close or unprocessed vertex
		FOR_CVF(f_it, trial_vh) {
			TriMesh::VertexHandle fixed_vh;
			TriMesh::VertexHandle close_or_unprocessed_vh;
			// Look for a fixed and a non-fixed vertex
			FOR_CFV(v_it, *f_it) {
				if(*v_it == trial_vh) continue;
				if(mesh.status(*v_it).tagged()) {
					fixed_vh = *v_it;
				} else {
					close_or_unprocessed_vh = *v_it;
				}
			}
			// (re-)calculcate distances from unprocessed and close vertices
			if(fixed_vh.is_valid() && close_or_unprocessed_vh.is_valid()) {
				try {
					double new_distance = compute_distance(mesh, close_or_unprocessed_vh, trial_vh, fixed_vh, geodesicDistanceProperty);
					update_distance(mesh, close_or_unprocessed_vh, new_distance, close_vertices, geodesicDistanceProperty);
				} catch(NoOverlap) {
					// Do not update the distance.
				}
			}
		}
	}
	for(auto vh : mesh.vertices()) {
		mesh.property(geodesicDistanceProperty, vh) = sqrt(mesh.property(geodesicDistanceProperty, vh));
	}
	mesh.release_vertex_status();
}
