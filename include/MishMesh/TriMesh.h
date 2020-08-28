#pragma once

#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <MishMesh/OpenMesh/DoublePrecisionTraits.h>

namespace MishMesh {
	// OpenMesh has since version 8.1 an own double precision traits implementation, but it
	// is incompatible to the DefaultTraits as it replaces the Vec3uc type for colors with
	// Vec3f. In addition, they still use Vec3f use for texture coordinates.
	typedef OpenMesh::TriMesh_ArrayKernelT<OpenMesh::DoublePrecisionTraits> TriMesh;
}
