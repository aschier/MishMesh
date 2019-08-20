#pragma once

namespace MishMesh{
	template<typename MeshT>
	void collapse_short_edges(MeshT &mesh, const double epsilon = 1e-8);
}
