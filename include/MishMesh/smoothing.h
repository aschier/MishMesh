#pragma once

#include <OpenMesh/Tools/Smoother/JacobiLaplaceSmootherT.hh>

namespace MishMesh {
	template<typename MeshT>
	using SmootherT = OpenMesh::Smoother::JacobiLaplaceSmootherT<MeshT>;
	template<typename MeshT>
	void smooth_mesh(MeshT &mesh, const int iterations = 1, const typename SmootherT<MeshT>::Component component = SmootherT<MeshT>::Component::Tangential_and_Normal);
}