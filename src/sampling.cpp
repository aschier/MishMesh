#include "MishMesh/sampling.h"

#include <MishMesh/TriMesh.h>
#include <MishMesh/utils.h>
#include <random>
#include <functional>
#include <cassert>

namespace MishMesh {
	static std::default_random_engine randgen(std::random_device{}());
	static auto random = [&](){return std::uniform_real_distribution<double>(0, 1)(randgen); };

	/**
	 * Sample a barycentric coordinate from a uniform random distribution.
	 * @param points Three points defining a triangle.
	 * @returns A barycentric coordinate.
	 */
	OpenMesh::Vec3d uniform_random_barycentric_coordinate(const std::array<OpenMesh::Vec3d, 3> points) {
		double xi1 = random();
		double xi2 = random();
		double s = sqrt(xi1);
		double t = xi2;
		auto bcoord = OpenMesh::Vec3d{1.0 - s, s*(1.0 - t), s*t};
		assert(bcoord.l1_norm() - 1.0 < FLT_EPSILON);
		return bcoord;
	}

	/**
	 * Sample a barycentric coordinate from a uniform random distribution.
	 * @param mesh A triangle mesh.
	 * @param fh A facehandle in the mesh.
	 * @returns The barycentric coordinate of the sampled point.
	 */
	OpenMesh::Vec3d uniform_random_barycentric_coordinate(const TriMesh & mesh, TriMesh::FaceHandle fh) {
		auto points = face_points(mesh, fh);
		return uniform_random_barycentric_coordinate(points);
	}

	/**
	 * Sample a random coordinate inside the triangle from a uniform random distribution.
	 * @param points Three points defining a triangle.
	 * @returns The coordinate of the sampled point.
	 */
	OpenMesh::Vec3d uniform_random_triangle_point(const std::array<OpenMesh::Vec3d, 3> points) {
		OpenMesh::Vec3d bcoord = uniform_random_barycentric_coordinate(points);
		return points[0] * bcoord[0] + points[1] * bcoord[1] + points[2] * bcoord[2];
	}

	/**
	 * Sample a random coordinate inside the triangle from a uniform random distribution.
	 * @param mesh A triangle mesh.
	 * @param fh A facehandle in the mesh.
	 * @returns The coordinate of the sampled point.
	 */
	OpenMesh::Vec3d uniform_random_triangle_point(const TriMesh &mesh, TriMesh::FaceHandle fh) {
		auto points = face_points(mesh, fh);
		return uniform_random_triangle_point(points);
	}

	/**
	 * Sample a random coordinate inside a annulus (ring) from a uniform random distribution.
	 * @param center The center of the annulus.
	 * @param inner_radius The inner radius.
	 * @param outer_radius The outer radius.
	 * @returns The coordinate of the sampled point.
	 */
	OpenMesh::Vec2d uniform_random_circle_point(const OpenMesh::Vec2d center, const double radius) {
		const double theta = 2 * M_PI * random();
		const double r = random() * radius;
		return {center[0] + r * cos(theta), center[1] + r * sin(theta)};
	}

	/**
	 * Sample a random coordinate inside a annulus (ring) from a uniform random distribution.
	 * @param center The center of the annulus.
	 * @param inner_radius The inner radius.
	 * @param outer_radius The outer radius.
	 * @returns The coordinate of the sampled point.
	 */
	OpenMesh::Vec2d uniform_random_annulus_point(const OpenMesh::Vec2d center, const double inner_radius, const double outer_radius) {
		const double theta = 2 * M_PI * random();
		const double r = sqrt(random() * (outer_radius*outer_radius - inner_radius * inner_radius) + inner_radius * inner_radius);
		return {center[0] + r * cos(theta), center[1] + r * sin(theta)};
	}
}
