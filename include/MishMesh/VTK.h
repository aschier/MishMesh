#pragma once

#include <MishMesh/TriMesh.h>

namespace MishMesh {
	typedef OpenMesh::VPropHandleT<double> VertexValueProperty;
	typedef OpenMesh::VPropHandleT<OpenMesh::Vec3d> VertexVectorProperty;
	typedef OpenMesh::FPropHandleT<double> FaceValueProperty;
	typedef OpenMesh::FPropHandleT<OpenMesh::Vec3d> FaceVectorProperty;

	void writeVTK(const MishMesh::TriMesh &mesh, const std::string filename,
	              VertexValueProperty prop_vertex_value = VertexValueProperty(),
	              VertexVectorProperty prop_vertex_vector = VertexVectorProperty(),
	              FaceValueProperty prop_face_value = FaceValueProperty(),
	              FaceVectorProperty prop_face_vector = FaceVectorProperty());
}