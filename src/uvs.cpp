#include "MishMesh/uvs.h"

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
			// Iterate around the to-vertex to find a mesh from the original polygon.
			auto heh2 = heh;
			do {
				heh2 = mesh.opposite_halfedge_handle(mesh.next_halfedge_handle(heh2));
				if(heh2.idx() < n_heh) {
					mesh.set_texcoord2D(heh, mesh.texcoord2D(heh2));
					break;
				}
			} while(heh != heh2);
		}
	}
}
