#pragma once

#include <MishMesh/TriMesh.h>
#include <OpenMesh/Tools/Smoother/JacobiLaplaceSmootherT.hh>

namespace MishMesh {
	typedef OpenMesh::Smoother::JacobiLaplaceSmootherT<TriMesh> SmootherT;
	void smooth_mesh(TriMesh &mesh, const int iterations = 1, const SmootherT::Component component = SmootherT::Component::Tangential_and_Normal);
}