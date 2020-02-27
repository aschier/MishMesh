#include "MishMesh/geodesics.h"
#include <MishMesh/macros.h>
#include <MishMesh/transformations.h>
#include <MishMesh/utils.h>
#include <OpenMesh/Tools/Utils/HeapT.hh>

using namespace std;
using namespace MishMesh;

inline OpenMesh::Vec2d rotate_ccw(const OpenMesh::Vec2d v) { return {-v[1], v[0]}; }
inline OpenMesh::Vec2d rotate_cw(const OpenMesh::Vec2d v) { return {v[1], -v[0]}; }

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
	assert(isfinite(edge_length_2));

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
	if(new_distance >= old_distance) return;

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

/*
 * Calculate section at an obtuse vertex, in which edges that split the triangle into two acute triangles can lie.
 * @param An array of 2D points defining the triangle.
 * @returns A pair of vectors defining the section.
 */
std::pair<OpenMesh::Vec2d, OpenMesh::Vec2d> calc_acute_section(const std::array<OpenMesh::Vec2d, 3> triangle_points) {
	assert(((triangle_points[0] - triangle_points[2]) | (triangle_points[1] - triangle_points[2])) < 0); // the angle is obtuse
	// Construct vectors orthogonal to the edges, that point inside the section spanned by
	// the triangle edges (v0, v2) and (v1, v2). The two vectors span an acute section, when the original section was obtuse.
	// We choose the indices such that vec0 belongs to the edge v0v2 and vec1 belongs to v1v2, because they point in the same direction.
	OpenMesh::Vec2d vec0 = rotate_cw(triangle_points[1] - triangle_points[2]);
	OpenMesh::Vec2d vec1 = rotate_ccw(triangle_points[0] - triangle_points[2]);

	// The section is acute
	assert((vec0 | vec1) > 0);

	// The boundary vectors point into the triangle
	assert(((triangle_points[0] - triangle_points[2]) | vec0) > 0);
	assert(((triangle_points[1] - triangle_points[2]) | vec1) > 0);

	return {vec0, vec1};
}

/**
 * Starting from an obtuse vertex, unfold triangles into the plane until a vertex is contained in the acute section of the obtuse vertex.
 * When this virtual vertex is found, the distance between the obtuse vertex and the virtual vertex can be calculated in the plane.
 * @param mesh The mesh.
 * @param heh The halfedge opposite of the obtuse vertex.
 * @param other_vh The third vertex in the triangle, that is not fixed, yet.
 * @param propDistance A GeodesicDistanceProperty for accessing the geodesic distances.
 * @returns a pair of the distance between the obtuse vertex and the virtual vertex in the plane and the
 *          vertex handle of the virtual vertex.
 */
std::pair<double, TriMesh::VertexHandle> find_virtual_vertex(TriMesh &mesh, TriMesh::HalfedgeHandle heh,
	const TriMesh::VertexHandle obtuse_vh, GeodesicDistanceProperty propDistance) {
	const auto points2D = embed_triangle(mesh, heh, mesh.point(obtuse_vh));
	auto trial_point_2D = points2D[0];
	const auto &obtuse_point_2D = points2D[2];

	auto trial_vh = mesh.from_vertex_handle(heh);
	auto other_vh = mesh.to_vertex_handle(heh);
	if(!mesh.status(trial_vh).tagged()) {
		std::swap(trial_vh, other_vh);
		trial_point_2D = points2D[1];
	}

	assert(obtuse_vertex(mesh, mesh.face_handle(heh)).is_valid()); // Vertex is obtuse in this triangle
	assert(mesh.status(trial_vh).tagged()); // trial vh is fixed
	assert(!mesh.status(other_vh).tagged()); // the other vertex is not fixed

	auto edge_vhs = edge_vertices(mesh, heh);

	// We assume an obtuse vertex v2, that is fixed and change the edge (v0, v1) by replacing
	// either v0 or v1 by v_opp in each iteration, until a fixed vertex v_opp is inside
	// the admissible sector. Then we use the edge (trial, v_opp) to calculate v2.
	//
	//           + v2 = obtuse
	//          / \
	//         /   \
	//        /     \
	// trial +-------+ other   // trial and other may be swapped
	//        \     / \
	//         \   /   \
	//          \ /     \     // Other unfolded triangles that already
	//           +-------+    // were processed in previous iterations
	//          / \     /
	//         /   \   /
	//        /     \ /
	//    v0 +-------+ v1
	//        \     /
	//         \   /
	//          \ /
	//           + v_opp

	// Calculate the sector that contains all possible bisectors, that split the triangle in two acute triangles
	OpenMesh::Vec2d t1, t2;
	std::tie(t1, t2) = calc_acute_section(points2D);
	// Normal vectors pointing inside the section.
	auto t1_orth = rotate_ccw(t1);
	auto t2_orth = rotate_cw(t2);

	// The points of the current virtual edge.
	auto edge_points = MishMesh::edge_points(mesh, heh);
	std::array<OpenMesh::Vec2d, 2> edge_points_2D{points2D[0], points2D[1]};
	// The new point below the current edge.
	OpenMesh::Vec2d opp_point_2D;

	MishMesh::TriMesh::VertexHandle virtual_vh = MishMesh::TriMesh::InvalidVertexHandle;
	while(!virtual_vh.is_valid()) {
		heh = mesh.opposite_halfedge_handle(heh);
		if(mesh.is_boundary(heh)) {
			return {numeric_limits<double>::infinity(), TriMesh::InvalidVertexHandle};
		}

		auto opp_vh = mesh.opposite_vh(heh); // the outgoing vertex "below" the edge
		auto opp_point_3D = mesh.point(opp_vh);

		// We do NOT use the vertices from the half edge heh, but the edge points from the previous iteration.
		// They are also not swapped, because we use opposite half edges only to find the next vertex
		// and not for orientation or vertex order.
		OpenMesh::Vec3d vh0_vh1 = edge_points[1] - edge_points[0]; // "horizontal" edge
		OpenMesh::Vec3d vh0_ovh = opp_point_3D - edge_points[0]; // new triangle vertex

		// Construct a orthogonal basis to project the new point into the 2D plane of the triangle.
		auto edge_vec2D_normalized = (edge_points_2D[1] - edge_points_2D[0]).normalized();
		OpenMesh::Vec2d edge_vec2D_orth_normalized = rotate_cw(edge_vec2D_normalized);

		assert(vh0_ovh.norm() > 0 && vh0_vh1.norm() > 0);
		double length_vh0_ovh_2 = vh0_ovh.sqrnorm();
		double l1 = vh0_vh1.normalized() | vh0_ovh; // parallel vector length
		double l1_2 = pow(vh0_vh1 | vh0_ovh, 2) / vh0_vh1.sqrnorm();
		double l2 = sqrt(length_vh0_ovh_2 - l1_2); // orthogonal vector length
		opp_point_2D = edge_points_2D[0] + l1 * edge_vec2D_normalized + l2 * edge_vec2D_orth_normalized;

		auto ivh_ovh_2D = opp_point_2D - obtuse_point_2D;

		if((ivh_ovh_2D | t1_orth) >= 0 && (ivh_ovh_2D | t2_orth) >= 0) {
			// The vertex is inside the angle.
			virtual_vh = opp_vh;
		} else if((ivh_ovh_2D | t2_orth) < 0) {
			// When the opp. vertex is right of the acute sector, we use the left edge
			// and the opp. vertex is the new v1
			heh = mesh.next_halfedge_handle(heh); // v0 -> v_opp
			edge_vhs[1] = opp_vh;
			edge_points[1] = opp_point_3D;
			edge_points_2D[1] = opp_point_2D;
		} else if((ivh_ovh_2D | t1_orth) < 0) {
			// When the opp. vertex is left of the acute sector, we use the right edge
			// and the opp. vertex is the new v0
			heh = mesh.prev_halfedge_handle(heh); // v_opp -> v1
			edge_vhs[0] = opp_vh;
			edge_points[0] = opp_point_3D;
			edge_points_2D[0] = opp_point_2D;
		} else {
			assert(false);
		}
	}

	if(!mesh.status(virtual_vh).tagged()) {
		// virtual vh is not fixed
		return {numeric_limits<double>::infinity(), TriMesh::InvalidVertexHandle};
	}

	double T1_2 = mesh.property(propDistance, trial_vh);
	assert(isfinite(T1_2));

	assert(obtuse_vh != edge_vhs[0] && obtuse_vh != edge_vhs[1]);
	double T2_2 = mesh.property(propDistance, virtual_vh);
	if(!isfinite(T1_2) || !isfinite(T2_2)) {
		return {numeric_limits<double>::infinity(), TriMesh::InvalidVertexHandle};
	}

	auto p1 = mesh.point(trial_vh);
	auto p2 = mesh.point(virtual_vh);

	try{
		// Origins in the reference coordinate system with edge (0,0), (length, 0)
		OpenMesh::Vec2d o1, o2;
		auto vec = (opp_point_2D - trial_point_2D);
		double vec_length_2 = vec.sqrnorm();
		double vec_length = vec.norm();
		//vec.normalize();
		OpenMesh::Vec2d vec_orth = rotate_ccw(vec);
		std::tie(o1, o2) = compute_projected_origins(vec_length_2, T1_2, T2_2);
		// Origins in the actual 2D coordinate system
		o1 = trial_point_2D + (o1[0] * vec + o1[1] * vec_orth) / vec_length;
		o2 = trial_point_2D + (o2[0] * vec + o2[1] * vec_orth) / vec_length;
		// Choose the nearer origin
		double T3_2 = max((points2D[2] - o1).sqrnorm(), (points2D[2] - o2).sqrnorm());
		assert(isfinite(T3_2));
		return {T3_2, virtual_vh};
	} catch(NoOverlap) {
		throw NoOverlap();
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
void MishMesh::compute_novotni_geodesics(TriMesh &mesh, const TriMesh::VertexHandle start_vh, const GeodesicDistanceProperty geodesicDistanceProperty, const bool handle_obtuse) {
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

	// Set the distance of the neighbor vertices from the edge length and change them from unprocessed to close
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
			TriMesh::VertexHandle close_or_unprocessed_vh1, close_or_unprocessed_vh2;
			// Look for a fixed and a non-fixed vertex
#ifdef _DEBUG
			int num_fixed = 0;
#endif
			FOR_CFV(v_it, *f_it) {
				if(*v_it == trial_vh) continue; // trial
				if(mesh.status(*v_it).tagged()) {
					fixed_vh = *v_it;
#ifdef _DEBUG
					num_fixed++;
#endif
					continue;
				}
				if(!close_or_unprocessed_vh1.is_valid()) {
					close_or_unprocessed_vh1 = *v_it;
				} else {
					assert(!close_or_unprocessed_vh2.is_valid());
					close_or_unprocessed_vh2 = *v_it;
				}
			}
			TriMesh::VertexHandle obtuse_vh = MishMesh::obtuse_vertex(mesh, *f_it);

			if(!close_or_unprocessed_vh1.is_valid()) {
				assert(!close_or_unprocessed_vh2.is_valid());
				// All vertices are fixed
				assert(num_fixed == 2);
				continue;
			}
			if(fixed_vh.is_valid()) {
				assert(!close_or_unprocessed_vh2.is_valid());
				assert(num_fixed == 1);
				// (re-)calculcate distances from unprocessed and close vertices
				try{
					double new_distance = compute_distance(mesh, close_or_unprocessed_vh1, trial_vh, fixed_vh, geodesicDistanceProperty);
					update_distance(mesh, close_or_unprocessed_vh1, new_distance, close_vertices, geodesicDistanceProperty);
				} catch(NoOverlap) {
					//continue;
				}
			} else if(handle_obtuse && obtuse_vh.is_valid() && (obtuse_vh == close_or_unprocessed_vh1 || obtuse_vh == close_or_unprocessed_vh2)) {
				assert(num_fixed == 0 || num_fixed == 1);
				// When the triangle contains an obtuse vertex that is not yet fixed,
				// calculate the corresponding virtual vertex and update the distance
				// calculated using the virtual vertex.
				assert(obtuse_vh != trial_vh && obtuse_vh != fixed_vh);

				try {
					TriMesh::VertexHandle other_vh = TriMesh::InvalidVertexHandle;
					FOR_CFV(v_it, *f_it) {
						if(*v_it != trial_vh && *v_it != obtuse_vh) {
							other_vh = *v_it;
							break;
						}
					}
					assert(other_vh.is_valid());
					assert(!mesh.status(other_vh).tagged());
					double virtual_distance;
					TriMesh::VertexHandle virtual_vh;
					auto heh = MishMesh::opposite_halfedge(mesh, *f_it, obtuse_vh);
					std::tie(virtual_distance, virtual_vh) = find_virtual_vertex(mesh, heh, obtuse_vh, geodesicDistanceProperty);
					update_distance(mesh, obtuse_vh, virtual_distance, close_vertices, geodesicDistanceProperty);
				} catch(NoOverlap) {
					// The vertices do not have the same projected origin
				}
			}
		}
	}
	for(auto vh : mesh.vertices()) {
		mesh.property(geodesicDistanceProperty, vh) = sqrt(mesh.property(geodesicDistanceProperty, vh));
	}
	mesh.release_vertex_status();
}
