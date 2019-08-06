#pragma once

#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <MishMesh/OpenMesh/DoublePrecisionTraits.h>

namespace MishMesh {
	typedef OpenMesh::PolyMesh_ArrayKernelT<OpenMesh::DoublePrecisionTraits> PolyMesh;
}
