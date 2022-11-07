/**
 * Create a dual sphere from a meshed sphere.
 * @tparam MeshT The mesh type, usually MishMesh::TriMesh or MishMesh::PolyMesh.
 * @param sphere_mesh A spherical mesh.
 * @note The face indices are expected to be continuous. You may need to call garbage_collect on the mesh
 *       before using this function.
*/
template<typename MeshT>
MishMesh::PolyMesh dualsphere(MeshT &sphere_mesh) {
    MishMesh::PolyMesh mesh;
    std::vector<typename MeshT::VertexHandle> dual_vertices(sphere_mesh.n_faces());

    for (auto fh : sphere_mesh.faces()) {
        std::vector<typename MeshT::VertexHandle> vhs;
        for(auto v_it=sphere_mesh.cfv_ccwbegin(fh); v_it != sphere_mesh.cfv_ccwend(fh); ++v_it) {
            vhs.push_back(*v_it);
        }
        auto p = (sphere_mesh.point(vhs[0]) + sphere_mesh.point(vhs[1]) + sphere_mesh.point(vhs[2])) / 3.0;
        p /= p.norm(); // normalize the sphere to radius 1
        dual_vertices[fh.idx()] = mesh.add_vertex(p);
    }
    for (auto vh : sphere_mesh.vertices()) {
        std::vector<MishMesh::PolyMesh::VertexHandle> vhs;
        for (auto f_it = sphere_mesh.cvf_ccwbegin(vh); f_it != sphere_mesh.cvf_ccwend(vh); f_it++) {
            vhs.push_back(dual_vertices[f_it->idx()]);
        }
        mesh.add_face(vhs.data(), vhs.size());
    }

    return mesh;
}
