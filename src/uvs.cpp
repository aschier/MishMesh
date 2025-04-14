#pragma once
#include "uvs.h"

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
}
