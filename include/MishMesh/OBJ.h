#pragma once

#include <MishMesh/TriMesh.h>
#include <MishMesh/PolyMesh.h>

#include <string>

namespace MishMesh {
	template<typename MeshT>
	bool writeOBJ(const MeshT &mesh, const std::string filename, int precision = 16);
	template<typename MeshT>
	bool readOBJ(MeshT &mesh, const std::string filename);
}
