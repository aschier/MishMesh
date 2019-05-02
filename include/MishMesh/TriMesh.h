#pragma once

#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <MishMesh/OpenMesh/DoublePrecisionTraits.h>

namespace MishMesh {
	typedef OpenMesh::TriMesh_ArrayKernelT<OpenMesh::DoublePrecisionTraits> TriMesh;
}
