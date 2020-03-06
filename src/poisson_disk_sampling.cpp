#include "MishMesh/poisson_disk_sampling.h"

#include <MishMesh/utils.h>
#include <MishMesh/BBox.h>
#include <MishMesh/sampling.h>

#include <deque>
#include <random>

namespace MishMesh {
	static std::default_random_engine randgen(std::random_device{}());
	inline int randint(int min, int max){
		return std::uniform_int_distribution<int>(min, max)(randgen);
	};

	/**
	 * A simple grid class, that allows to use two-dimensional indices with a std::vector as storage class.
	 * @tparam T the data type to store.
	 */
	template<typename T>
	class BackgroundGrid{
	public:
		BackgroundGrid(uint width, uint height, T initial): width(width), height(height), data(std::vector<T>(width*height, initial)) {};
		inline void set(int x, int y, const T &value) { data[width*y + x] = value; }
		inline const T operator()(int x, int y) const { return data[width*y + x]; }
	private:
		uint width, height;
		std::vector<T> data;
	};

	/**
	 * Find the grid indices (i, j) of a point using the grid dimensions and the number of grid cells.
	 * @param bbox A bounding box, that contains the grid
	 * @param Nx the number of grid cells in x direction.
	 * @param Ny the number of grid cells in y direction.
	 * @param point The input point.
	 * @returns A pair of indices of grid cells.
	 */
	inline std::pair<uint, uint> grid_idx(const MishMesh::BBox<OpenMesh::Vec2d, 2> &bbox, uint Nx, uint Ny, const OpenMesh::Vec2d &point){
		assert(bbox.contains(point));
		const double width = bbox.rbn[0] - bbox.ltf[0];
		const double height = bbox.rbn[1] - bbox.ltf[1];
		std::pair<uint, uint> result;
		result = std::make_pair(
			static_cast<uint>(((point[0] - bbox.ltf[0])) / width * Nx),
			static_cast<uint>(((point[1] - bbox.ltf[1])) / height * Ny)
		);
		assert(result.first <= Nx && result.second <= Ny);
		return result;
	}

	/**
	 * Test if a point is in a triangle given by three points.
	 * @param triangle_points Three points, that define a triangle.
	 * @param p The point
	 * @returns true, when the point lies inside the triangle.
	 */
	inline bool point_in_triangle(const std::array<OpenMesh::Vec2d, 3> &triangle_points, const OpenMesh::Vec2d &p) {
		OpenMesh::Vec2d BA = triangle_points[1] - triangle_points[0];
		OpenMesh::Vec2d CA = triangle_points[2] - triangle_points[0];
		OpenMesh::Vec2d PA = p - triangle_points[0];
		double d00 = CA | CA;
		double d01 = CA | BA;
		double d02 = CA | PA;
		double d11 = BA | BA;
		double d12 = BA | PA;

		double denominator = d00 * d11 - d01 * d01;
		double u = (d11 * d02 - d01 * d12) / denominator;
		double v = (d00 * d12 - d01 * d02) / denominator;
		return u > 0 && v > 0 && (u + v) < 1.0;
	}

	/**
	 * An internal function for a point in triangle test, when certain values used in the algorithm
	 * are already calculated.
	 * @note Do not use this function directly, but use get_fast_point_in_triangle.
	 * @see get_fast_point_in_triangle_function
	 */
	inline bool _fast_point_in_triangle(
		const OpenMesh::Vec2d p0, const OpenMesh::Vec2d BA, const OpenMesh::Vec2d CA,
		const double d00, const double d01, const double d11, const double denominator,
		const OpenMesh::Vec2d &p) {
		OpenMesh::Vec2d PA = p - p0;
		double d02 = CA | PA;
		double d12 = BA | PA;

		double u = (d11 * d02 - d01 * d12) / denominator;
		double v = (d00 * d12 - d01 * d02) / denominator;
		return u > 0 && v > 0 && (u + v) < 1.0;
	}

	/**
	 * Get a function that computes if a given point lies inside the triangle used to create the function.
	 * This function precomputes some values that only depend on the triangle itself and not on the input point,
	 * such that the returned function does not need to compute them in each call.
	 * @param triangle_points Three points, that define a triangle.
	 * @returns A function bool point_in_triangle(OpenMesh::Vec2d point), that tests if the point lies inside the
	 *          triangle given by triangle_points.
	 */
	std::function<bool(const OpenMesh::Vec2d &p)> get_fast_point_in_triangle_function(const std::array<OpenMesh::Vec2d, 3> &triangle_points) {
		OpenMesh::Vec2d BA = triangle_points[1] - triangle_points[0];
		OpenMesh::Vec2d CA = triangle_points[2] - triangle_points[0];
		double d00 = CA | CA;
		double d01 = CA | BA;
		double d11 = BA | BA;

		double denominator = d00 * d11 - d01 * d01;
		return std::bind(_fast_point_in_triangle, triangle_points[0], BA, CA, d00, d01, d11, denominator, std::placeholders::_1);
	}

	/**
	 * Subsample points in a circle using poisson disk sampling using the algorithm of Bridson.
	 *
	 * @param center The center of the circle.
	 * @param radius The radius of the circle.
	 * @param min_distance The minimum distance between two points.
	 * @param cell_size The cell size for a background grid, that is used for checking if there is a point
	 *                  with distance less than min_distance.
	 * @param max_tries The number of tries to find a new point in the neighborhood of an active point,
	 *                  before continuing with next active point.
	 * @returns a vector containing the sampled points.
	 * @note For the algorithm, you need to make sure that cell_size <= min_distance / sqrt(2)
	 */
	std::vector<OpenMesh::Vec2d> poisson_disk_sampling_circle(const OpenMesh::Vec2d center, const double radius, const double min_distance, const double cell_size, const int max_tries) {
		BBox<OpenMesh::Vec2d, 2> bbox({center[0] - radius, center[1] - radius}, {center[0] + radius, center[1] + radius});
		const double radius_2 = radius * radius;
		std::function<bool(const OpenMesh::Vec2d &)> point_in_sphere_test = [&](const OpenMesh::Vec2d &p) { return (center - p).sqrnorm() < radius_2; };
		return poisson_disk_sampling(bbox, point_in_sphere_test, {center}, min_distance, cell_size, max_tries);
	}

	/**
	 * Subsample points in a triangle using poisson disk sampling using the algorithm of Bridson.
	 *
	 * @param triangle_points Three points, that define a triangle.
	 * @param min_distance The minimum distance between two points.
	 * @param cell_size The cell size for a background grid, that is used for checking if there is a point
	 *                  with distance less than min_distance.
	 * @param max_tries The number of tries to find a new point in the neighborhood of an active point,
	 *                  before continuing with next active point.
	 * @returns a vector containing the sampled points.
	 * @note For the algorithm, you need to make sure that cell_size <= min_distance / sqrt(2)
	 */
	std::vector<OpenMesh::Vec2d> poisson_disk_sampling(const std::array<OpenMesh::Vec2d, 3> triangle_points, const double min_distance, const double cell_size, const int max_tries) {
		auto bbox = bounding_box<2, OpenMesh::Vec2d>(triangle_points.begin(), triangle_points.end());
		auto fast_point_in_triangle = get_fast_point_in_triangle_function(triangle_points);
		return poisson_disk_sampling(bbox, fast_point_in_triangle, std::vector<OpenMesh::Vec2d>(triangle_points.begin(), triangle_points.end()), min_distance, cell_size, max_tries);
	}

	/**
	 * Subsample points using poisson disk sampling using the algorithm of Bridson.
	 * [Bridson, R. (2007). Fast Poisson disk sampling in arbitrary dimensions.
	 * In ACM SIGGRAPH 2007 Sketches on - SIGGRAPH ’07, (San Diego, California: ACM Press), pp. 22-es.]
	 * https://www.cs.ubc.ca/~rbridson/docs/bridson-siggraph07-poissondisk.pdf
	 *
	 * @param bbox The bounding box of the target shape.
	 * @param point_test A function, that returns true, when the point lies inside the target shape.
	 * @param initial_points Points that are already given before sampling.
	 * @param min_distance The minimum distance between two points.
	 * @param cell_size The cell size for a background grid, that is used for checking if there is a point
	 *                  with distance less than min_distance.
	 * @param max_tries The number of tries to find a new point in the neighborhood of an active point,
	 *                  before continuing with next active point.
	 * @returns a vector containing the sampled points.
	 * @note For the algorithm, you need to make sure that cell_size <= min_distance / sqrt(2)
	 */
	std::vector<OpenMesh::Vec2d> poisson_disk_sampling(MishMesh::BBox<OpenMesh::Vec2d, 2> bbox, std::function<bool(const OpenMesh::Vec2d &)> point_test, const std::vector<OpenMesh::Vec2d> initial_points,
		const double min_distance, const double cell_size, const int max_tries) {
		std::vector<OpenMesh::Vec2d> result;
		assert(cell_size <= min_distance / sqrt(2));
		const double max_length = bbox.max_side_length();
		const uint N = static_cast<uint>(max_length / cell_size);
		BackgroundGrid<bool> grid(N + 1, N + 1, false);

		auto _grid_idx = std::bind(grid_idx, bbox, N, N, std::placeholders::_1);
		uint idx_x, idx_y;

		std::vector<OpenMesh::Vec2d> active_list;
		for(auto &p : initial_points){
			std::tie(idx_x, idx_y) = _grid_idx(p);
			grid.set(idx_x, idx_y, 1);
			active_list.push_back(p);
		}

		while(!active_list.empty()) {
			// Choose a random point and remove it from the list
			int idx = randint(0, static_cast<int>(active_list.size() - 1));
			auto p = active_list[idx];
			active_list[idx] = active_list.back();
			active_list.pop_back();
			for(int t = 0; t < max_tries; t++) {
				auto new_point = MishMesh::uniform_random_annulus_point(p, min_distance * 2, min_distance);
				if(!point_test(new_point)) continue;
				std::tie(idx_x, idx_y) = _grid_idx(new_point, true);
				bool admissible = true;

				auto from_point = bbox.clip(new_point - OpenMesh::Vec2d{min_distance, min_distance});
				auto grid_from = _grid_idx(from_point);
				auto to_point = bbox.clip(new_point + OpenMesh::Vec2d{min_distance, min_distance});
				auto grid_to = _grid_idx(to_point);
				// When clipping at the boundary, both points may lie in the same (boundary) cell.
				if(grid_from == grid_to) continue;

				// When one of the cells in the grid around the new point is occupied by another point,
				// the new point is not admissible.
				for(uint i = grid_from.first; i < grid_to.first; i++) for(uint j = grid_from.second; j < grid_to.second; j++) {
					if(grid(i, j)){
						admissible = false;
						break;
					}
				}
				if(!admissible) continue;

				// Store the point and mark it as occupied in the background grid.
				grid.set(idx_x, idx_y, true);
				active_list.push_back(new_point);
				result.push_back(new_point);
			}
		}
		return result;
	}

	/**
	 * Subsample points in a triangle using poisson disk sampling.
	 * @param triangle_points Three points, that define a triangle.
	 * @param min_distance The minimum distance between two points.
	 * @returns a vector containing the sampled points.
	 */
	std::vector<OpenMesh::Vec2d> poisson_disk_sampling(const std::array<OpenMesh::Vec2d, 3> triangle_points, const double min_distance){
		return poisson_disk_sampling(triangle_points, min_distance, min_distance / sqrt(2.0));
	};
}
