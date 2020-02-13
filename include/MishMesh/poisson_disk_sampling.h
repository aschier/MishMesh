#pragma once

#include <MishMesh/TriMesh.h>
#include <MishMesh/BBox.h>
#include <functional>

namespace MishMesh {
	std::vector<OpenMesh::Vec2d> poisson_disk_sampling_circle(const OpenMesh::Vec2d center, const double radius, const double min_distance, const double cell_size, const int max_tries = 30);
	std::vector<OpenMesh::Vec2d> poisson_disk_sampling(const std::array<OpenMesh::Vec2d, 3> triangle_points, const double min_distance, const double cell_size, const int max_tries = 30);
	std::vector<OpenMesh::Vec2d> poisson_disk_sampling(MishMesh::BBox<OpenMesh::Vec2d, 2> bbox, std::function<bool(const OpenMesh::Vec2d&)> point_test, const std::vector<OpenMesh::Vec2d> initial_points, const double min_distance, const double cell_size, const int max_tries);
	std::vector<OpenMesh::Vec2d> poisson_disk_sampling(const std::array<OpenMesh::Vec2d,3> triangle_points, const double min_distance);
}