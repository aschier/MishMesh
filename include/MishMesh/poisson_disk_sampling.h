#pragma once

#include <MishMesh/TriMesh.h>

namespace MishMesh {
	std::vector<OpenMesh::Vec2d> poisson_disk_sampling(const std::array<OpenMesh::Vec2d, 3> triangle_points, const double min_distance, const double cell_size, const int max_tries = 30);
	std::vector<OpenMesh::Vec2d> poisson_disk_sampling(const std::array<OpenMesh::Vec2d,3> triangle_points, const double min_distance);
}