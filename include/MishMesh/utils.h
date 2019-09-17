#pragma once
#include <MishMesh/TriMesh.h>
#include <array>
#include <limits>

namespace MishMesh {
	TriMesh::HalfedgeHandle opposite_halfedge(const TriMesh &mesh, const TriMesh::FaceHandle &fh, TriMesh::VertexHandle &vh);
	std::array<TriMesh::VertexHandle, 3> face_vertices(const TriMesh &mesh, const TriMesh::FaceHandle fh);

	template<int DIM>
	double compute_area(const std::array<OpenMesh::VectorT<double, DIM>, 3> points);
	double compute_area(const TriMesh &mesh, const std::array<TriMesh::VertexHandle, 3> vertices);
	double compute_area(const TriMesh &mesh, const TriMesh::FaceHandle fh);

	/** A simple bounding box defined by the left-top(-far) and right-bottom(-near) point
	 * @tparam VectorT A vector class, that stores DIM entries and provides an operator[] for the entries.
	 * @tparam DIM The dimension of the point and the bounding box.
	 */
	template<typename VectorT = OpenMesh::Vec3d, unsigned int DIM = 3>
	struct BBox {
		VectorT ltf;
		VectorT rbn;
		BBox(){
			for(short j = 0; j < DIM; j++){
				ltf[j] = 0.0;
				rbn[j] = 0.0;
			}
		}

		/**
		 * Create a BBox with all coordinates 0.
		 */
		static BBox zero() {
			BBox bbox;
			for(short j = 0; j < DIM; j++){
				bbox.ltf[j] = 0.0;
				bbox.rbn[j] = 0.0;
			}
			return bbox;
		}

		/**
		 * Create a BBox between -infinity and +infinity in all dimensions.
		 */
		static BBox infinity() {
			BBox bbox;
			for(short j = 0; j < DIM; j++){
				bbox.ltf[j] = std::numeric_limits<double>::infinity();
				bbox.rbn[j] = -std::numeric_limits<double>::infinity();
			}
			return bbox;
		}

		BBox(VectorT ltf, VectorT rbn): ltf(ltf), rbn(rbn) {};
		/**
		 * A BBox is valid, if all left-top(-far) values are smaller than the corresponding right-bottom(-near) values.
		 * @returns true, when the bounding box is valid.
		 */
		bool is_valid() const {
			for(unsigned int j = 0; j < DIM; j++){
				if(ltf[j] > rbn[j]) {
					return false;
				}
			}
			return true;
		}
		/**
		 * Make the bounding box valid, by swapping min/max coordinates, that are in the wrong order,
		 * i.e. when bbox.ltf[j] < bbox.rbn[j], bbox.ltf[j] and bbox.rbn[j] are swapped.
		 */
		void make_valid() {
			for(unsigned int j = 0; j < DIM; j++){
				if(ltf[j] > rbn[j]) {
					std::swap(ltf[j], rbn[j]);
				}
			}
		}
		/**
		 * Check if a point is in a given bounding box.
		 * @param point the point.
		 * @param bbox The bounding box.
		 * @param exact When exact is true, the point must be inside the bounding box,
		 *        otherwise it must be inside a bounding box that is FLT_EPSILON in each direction.
		 * @returns true, when the point is inside the bounding box and false otherwise.
		 */
		bool contains(const VectorT point, bool exact = false) const {
			double padding = exact ? 0 : FLT_EPSILON;
			for(unsigned int j = 0; j < DIM; j++){
				if(point[j] > rbn[j] + padding || point[j] < ltf[j] - padding){
					return false;
				}
			}
			return true;
		}
		/**
		 * Clip a point to a bounding box.
		 * @param point The point.
		 * @returns The point, if it is inside the bounding box or a clipped point, when the point was outside the box.
		 */
		VectorT clip(VectorT point) const {
			for(unsigned int j = 0; j < DIM; j++){
				if(point[j] > rbn[j]){
					point[j] = rbn[j];
				} else if(point[j] < ltf[j]){
					point[j] = ltf[j];
				}
			}
			return point;
		}

		double diameter() const {
			double result = 0.0;
			for(short j = 0; j < DIM; j++) {
				result += (rbn[j] - ltf[j])*(rbn[j] - ltf[j]);
			}
			return sqrt(result);
		}
	};

	/**
	 * Calculate the axis aligned bounding box of a mesh.
	 * @param mesh The mesh.
	 * @tparam MeshT The mesh type.
	 * @tparam A bbox type with a vector that has 3 components.
	 */
	template<typename MeshT, typename BBoxT = MishMesh::BBox<OpenMesh::Vec3d, 3>>
	BBoxT bounding_box(const MeshT &mesh) {
		BBoxT result;
		for(short j = 0; j < 3; j++) {
			result.ltf[j] = std::numeric_limits<double>::infinity();
			result.rbn[j] = -std::numeric_limits<double>::infinity();
		}
		for(auto vh : mesh.vertices()) {
			auto &p = mesh.point(vh);
			for(short j = 0; j < 3; j++) {
				result.ltf[j] = std::min(result.ltf[j], p[j]);
				result.rbn[j] = std::max(result.rbn[j], p[j]);
			}
		}
		return result;
	}
}
