#include <OpenMesh/Core/Geometry/VectorT.hh>
#include <MishMesh/BBox.h>

namespace MishMesh {
	namespace Meshing {
		enum FaceDirection {
			LEFT,
			TOP,
			FAR
		};

		class Grid {
		public:
			virtual int resolution(int dim) const = 0;
			virtual OpenMesh::Vec3i resolution() const = 0;
			virtual OpenMesh::Vec3i numCells() const = 0;
			virtual OpenMesh::Vec3d point(int i, int j, int k) const = 0;
			virtual OpenMesh::Vec3d point(int idx) const = 0;
			virtual OpenMesh::Vec3d point(OpenMesh::Vec3i idxs) const = 0;
			/// Get a point index from coordinates
			virtual int point_idx(int i, int j, int k) const = 0;
			/// Get a point index relative to another point index, i.e.,
			/// i steps in x, j steps in y and k steps in z direction
			virtual int point_idx(int current_idx, int i, int j, int k) const = 0;
			virtual OpenMesh::Vec3i grid_idxs(int idx) const = 0;

			/**
			 * Generate the grid indices for a face given by a grid point (left, top, far) and
			 * a dimension
			 * @param faceDirection The direction of the cube face (left, top, far).
			 * @param i x coordinate.
			 * @param j y coordinate.
			 * @param k z coordinate.
			 */
			inline std::array<int, 4> face_idxs(const FaceDirection faceDirection, const int i, const int j, const int k) const {
				if(faceDirection == FaceDirection::LEFT) {
					// left face
					return std::array<int, 4>{
					    point_idx(i, j, k),
					    point_idx(i, j, k + 1),
					    point_idx(i, j + 1, k + 1),
					    point_idx(i, j + 1, k),
					};
				} else if(faceDirection == FaceDirection::TOP) {
					// top face
					return std::array<int, 4>{
					    point_idx(i, j, k),
					    point_idx(i + 1, j, k),
					    point_idx(i + 1, j, k + 1),
					    point_idx(i, j, k + 1),
					};
				} else {
					// far face
					return std::array<int, 4>{
					    point_idx(i, j, k),
					    point_idx(i, j + 1, k),
					    point_idx(i + 1, j + 1, k),
					    point_idx(i + 1, j, k),
					};
				}
			}
		};

		class RegularGrid : public Grid {
		public:
			inline RegularGrid(){};
			inline RegularGrid(MishMesh::BBox<OpenMesh::Vec3d, 3> bbox, OpenMesh::Vec3i steps) : bbox(bbox), steps(steps){};
			inline int resolution(int dim) const override {
				return steps[dim];
			};
			inline OpenMesh::Vec3i resolution() const override {
				return steps;
			};
			inline OpenMesh::Vec3i numCells() const override {
				return resolution() - OpenMesh::Vec3i{1, 1, 1};
			};
			virtual OpenMesh::Vec3d point(int i, int j, int k) const override {
				return OpenMesh::Vec3d{
				    bbox.ltf[0] + (bbox.rbn[0] - bbox.ltf[0]) * i / (steps[0] - 1),
				    bbox.ltf[1] + (bbox.rbn[1] - bbox.ltf[1]) * j / (steps[1] - 1),
				    bbox.ltf[2] + (bbox.rbn[2] - bbox.ltf[2]) * k / (steps[2] - 1)};
			};

			/**
				 *  Get the index of a grid point.
				 * @param i The x index.
				 * @param j The y index.
				 * @param k The z index.
				 * @returns The index of the grid point.
				 */
			inline int point_idx(int i, int j, int k) const override {
				return i + resolution(0) * j + resolution(0) * resolution(1) * k;
			}
			inline OpenMesh::Vec3d point(OpenMesh::Vec3i idxs) const override {
				return point(idxs[0], idxs[1], idxs[2]);
			};
			inline OpenMesh::Vec3d point(int idx) const override {
				return point(grid_idxs(idx));
			};

			/**
				 * Get the index of a grid point relative to a given index.
				 * @param current_idx The index of a grid point.
				 * @param i The x offset relative to the given point.
				 * @param j The y offset relative to the given point.
				 * @param k The z offset relative to the given point.
				 * @returns The index of the new point.
				 */
			inline int point_idx(int current_idx, int i, int j, int k) const override {
				return current_idx + point_idx(i, j, k);
			}

			/**
				 * Generate the grid indices for given point index.
				 * @param idx The point index.
				 * @returns The x, y, z indices in the grid.
				 */
			inline OpenMesh::Vec3i grid_idxs(int idx) const override {
				const int res_x = resolution(0);
				const int res_y = resolution(1);
				const int res_z = resolution(2);
				return OpenMesh::Vec3i{
				    idx % res_x,
				    idx / res_x % res_y,
				    idx / (res_x * res_y)};
			}

		private:
			OpenMesh::Vec3i steps;
			BBox<OpenMesh::Vec3d, 3> bbox;
		};
	}
}
