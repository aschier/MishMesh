#include <MishMesh/utils.h>
#include <MishMesh/macros.h>

using namespace std;

namespace MishMesh {
	TriMesh::HalfedgeHandle opposite_halfedge(const TriMesh &mesh, const TriMesh::FaceHandle &fh, TriMesh::VertexHandle &vh) {
		TriMesh::HalfedgeHandle result_heh;
		FOR_CFH(h_it, fh) {
			if(mesh.from_vertex_handle(*h_it) != vh && mesh.to_vertex_handle(*h_it) != vh) {
				result_heh = *h_it;
				break;
			}
		}
		assert(result_heh.is_valid());
		return result_heh;
	}

	std::array<TriMesh::VertexHandle, 3> face_vertices(const TriMesh &mesh, const TriMesh::FaceHandle fh) {
		array<TriMesh::VertexHandle, 3> vhs;
		short j = 0;
		FOR_CFV(v_it, fh) {
			vhs[j++] = *v_it;
		}
		return vhs;
	}
}
