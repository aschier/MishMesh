#include "MishMesh/uvs.h"
#include "MishMesh/utils.h"

namespace MishMesh {
	/**
	 * Flip an edge in the mesh and adapt the UV coordinates of the halfedges to the
	 * UV coordinates at the new vertices.
	 *
	 * @param mesh The mesh.
	 * @param eh The edge handle of the edge to be flipped.
	 * @note The function assumes that the mesh has halfedge texture coordinates
	 *       and that eh is no UV boundary.
	 */
	void flip_with_uvs(MishMesh::TriMesh &mesh, const MishMesh::TriMesh::EdgeHandle &eh) {
		mesh.flip(eh);
		auto heh1 = mesh.halfedge_handle(eh, 0);
		auto heh2 = mesh.opposite_halfedge_handle(heh1);
		mesh.set_texcoord2D(heh1, mesh.texcoord2D(mesh.prev_halfedge_handle(mesh.opposite_halfedge_handle(heh1))));
		mesh.set_texcoord2D(heh2, mesh.texcoord2D(mesh.prev_halfedge_handle(mesh.opposite_halfedge_handle(heh2))));
	}

	/**
	 * Split an edge and interpolate the UV coordinates for all halfedges adjacent to the split
	 * vertex. The UVs at the halfedge and its opposite halfedge are interpolated separately,
	 * such that UV boundaries are preserved.
	 *
	 * The edge and halfedge properties are copied to the new segments using mesh.copy_all_properties
	 * and for the inserted triangulation edges and new faces all properties but the UV coordinates are
	 * undefined.
	 *
	 * @param mesh The mesh.
	 * @param heh The boundary halfedge handle to split.
	 * @param t the relative position in [0, 1] where to split the edge.
	 * @returns A pair of halfedge handles representing the two segments of the input halfedge.
	 * @note the mesh is assumed to have a texcoord2D halfedge property that stores the
	 *       current UVs.
	 * @note The original halfedge becomes one of the two segments. The new vertex can be retrieved
	 *       with result.first.to().
	 */
	SmartHalfedgePair split_with_uvs(MishMesh::TriMesh &mesh, const OpenMesh::SmartHalfedgeHandle heh, const double t) {
		const auto p = (1.0 - t) * mesh.point(heh.from()) + t * mesh.point(heh.to());
		auto vh = mesh.add_vertex(p);

		auto interpolate_uv = [&](const auto &heh, const double t) {
			return (1.0 - t) * mesh.texcoord2D(heh.prev()) + t * mesh.texcoord2D(heh);
		};

		MishMesh::TriMesh::TexCoord2D edge_from_uv = mesh.texcoord2D(heh.opp());
		MishMesh::TriMesh::TexCoord2D edge_to_uv = mesh.texcoord2D(heh);
		MishMesh::TriMesh::TexCoord2D path_from_uv = mesh.texcoord2D(heh.next());
		MishMesh::TriMesh::TexCoord2D path_to_uv = mesh.texcoord2D(heh.opp().next());
		MishMesh::TriMesh::TexCoord2D interpolated_uv = interpolate_uv(heh, t);
		MishMesh::TriMesh::TexCoord2D interpolated_uv_opp = interpolate_uv(heh.opp(), 1 - t);

		auto [edge_heh1, edge_heh2] = MishMesh::split_edge(mesh, heh, vh);
		const auto &new_heh = (edge_heh1 == heh) ? edge_heh2 : edge_heh1;

		mesh.copy_all_properties(heh.edge(), new_heh.edge());
		mesh.copy_all_properties(heh, new_heh);
		mesh.copy_all_properties(heh.opp(), new_heh.opp());

		// Direction: current_vh = edge_heh2.prev().from() -> new_vh -> edge_heh2.opp().next().to()
		auto triangulation_heh1 = edge_heh2.prev();
		auto triangulation_heh2 = edge_heh2.opp().next();

		mesh.set_texcoord2D(edge_heh1, interpolated_uv);
		mesh.set_texcoord2D(edge_heh1.opp(), edge_from_uv);
		mesh.set_texcoord2D(edge_heh2, edge_to_uv);
		mesh.set_texcoord2D(edge_heh2.opp(), interpolated_uv_opp);
		if(!edge_heh2.is_boundary()) {
			mesh.set_texcoord2D(triangulation_heh1, interpolated_uv);
			mesh.set_texcoord2D(triangulation_heh1.opp(), path_from_uv);
		}
		if(!edge_heh2.opp().is_boundary()) {
			mesh.set_texcoord2D(triangulation_heh2, path_to_uv);
			mesh.set_texcoord2D(triangulation_heh2.opp(), interpolated_uv_opp);
		}
		return {edge_heh1, edge_heh2};
	}

	/**
	 * Triangulate a PolyMesh and assign the new halfedges inside the faces the UV coordinates of
	 * the boundary halfedges of the original polygonal faces.
	 *
	 * @param[inout] mesh The mesh.
	 * @note The mesh is triangulated in-place using mesh.triangulate().
	 */
	void triangulate_with_uvs(MishMesh::PolyMesh &mesh) {
		int n_heh = mesh.n_halfedges(); // PolyMesh halfedges
		mesh.triangulate();
		// New halfedges
		for(int i = n_heh; i < mesh.n_halfedges(); i++) {
			MishMesh::PolyMesh::HalfedgeHandle heh = mesh.halfedge_handle(i);
			// Iterate around the to-vertex to find an edge from the original polygon.
			auto heh2 = heh;
			do {
				heh2 = mesh.prev_halfedge_handle(mesh.opposite_halfedge_handle(heh2));
				if(heh2.idx() < n_heh) {
					mesh.set_texcoord2D(heh, mesh.texcoord2D(heh2));
					break;
				}
			} while(heh != heh2);
		}
	}
}
