#pragma once

#include <MishMesh/TriMesh.h>
#include <MishMesh/PolyMesh.h>

namespace MishMesh {
	typedef OpenMesh::VPropHandleT<double> VertexValueProperty;
	typedef OpenMesh::VPropHandleT<OpenMesh::Vec3d> VertexVectorProperty;
	typedef OpenMesh::FPropHandleT<double> FaceValueProperty;
	typedef OpenMesh::FPropHandleT<OpenMesh::Vec3d> FaceVectorProperty;

	template<typename MeshT>
	void writeVTK(const MeshT &mesh, const std::string filename, bool is_trimesh,
	                        VertexValueProperty prop_vertex_value,
	                        VertexVectorProperty prop_vertex_vector,
	                        FaceValueProperty prop_face_value,
	                        FaceVectorProperty prop_face_vector);

	void writeVTK(const MishMesh::TriMesh &mesh, const std::string filename,
	              VertexValueProperty prop_vertex_value = VertexValueProperty(),
	              VertexVectorProperty prop_vertex_vector = VertexVectorProperty(),
	              FaceValueProperty prop_face_value = FaceValueProperty(),
	              FaceVectorProperty prop_face_vector = FaceVectorProperty());

	void writeVTK(const MishMesh::PolyMesh &mesh, const std::string filename,
	              VertexValueProperty prop_vertex_value = VertexValueProperty(),
	              VertexVectorProperty prop_vertex_vector = VertexVectorProperty(),
	              FaceValueProperty prop_face_value = FaceValueProperty(),
	              FaceVectorProperty prop_face_vector = FaceVectorProperty());
}