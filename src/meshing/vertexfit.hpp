#include <MishMesh/PolyMesh.h>

#ifdef HAS_EIGEN
#include <Eigen/Eigen>
#endif

/*
 * Fit vertices to a plane defined by a set of points with normals as described in
 * Ju, Tao, et al. "Dual contouring of hermite data." Proceedings of the 29th annual conference on Computer graphics and interactive techniques. 2002.
 * and
 * Lindstrom, P. (2000, July). Out-of-core simplification of large polygonal models. In Proceedings of the 27th annual conference on Computer graphics and interactive techniques (pp. 259-262).
 */

namespace MishMesh {
	namespace Meshing {
		enum SVDCutoffType {
			RELATIVE,
			ABSOLUTE
		};
		constexpr double SVD_CUTOFF_FACTOR = 0.1;

		/**
		 * Fit a point to the average of a set of points and clip it to a bounding box.
		 * @param The mesh.
		 * @param vh The vertex handle for the point to fit.
		 * @param points a set of points, that is used to calculcate the average.
		 */
		void fit_vertex_average(MishMesh::PolyMesh &mesh,
		                        const MishMesh::PolyMesh::VertexHandle vh,
		                        const MishMesh::BBox<OpenMesh::Vec3d, 3> &bbox,
		                        const std::vector<OpenMesh::Vec3d> &points) {
			OpenMesh::Vec3d x = std::accumulate(points.begin(), points.end(), OpenMesh::Vec3d{0, 0, 0}) / points.size();
			mesh.set_point(vh, OpenMesh::Vec3d(x.data()));
		}

#ifdef HAS_EIGEN
		/**
		 * Fit a point to a plane defined by a set of points with normals.
		 * @param The mesh.
		 * @param vh The vertex handle for the point to fit.
		 * @param bbox a bounding box that is used to clip the point when it lies outside.
		 * @param points a set of points, that is used to calculcate the average.
		 * @param normals The normals at the points in the points array.
		 * @param clip If clip is true, the result will be clipped on the bounding box.
		 * @param svd_cutoff_type For absolute cutoff, singular values smaller than svd_cutoff_factor
		 *        are truncated and for relative cutoff singular values S_i with S_i/S_0 smaller than
		 *        svd_cutoff_factor are truncated.
		 * @param svd_cutoff_factor Threshold which singular values are truncated to avoid instable results,
		 *        when determining the plane of the points.
		 */
		void fit_vertex_svd(MishMesh::PolyMesh &mesh,
		                    const MishMesh::PolyMesh::VertexHandle vh,
		                    const MishMesh::BBox<OpenMesh::Vec3d, 3> &bbox,
		                    const std::vector<OpenMesh::Vec3d> &points,
		                    const std::vector<OpenMesh::Vec3d> &normals,
		                    const bool clip,
		                    const SVDCutoffType svd_cutoff_type,
		                    const double svd_cutoff_factor) {
			Eigen::Matrix3d ATA;
			Eigen::Vector3d ATb;
			ATA.setZero();
			ATb.setZero();
			assert(points.size() == normals.size());
			MishMesh::BBox<Eigen::Vector3d, 3> cube_bbox_eigen{Eigen::Vector3d(bbox.ltf.data()), Eigen::Vector3d(bbox.rbn.data())};
			for(int i = 0; i < points.size(); i++) {
				Eigen::Vector3d p_i(points[i].data());
				Eigen::Matrix<double, 3, 1> n_i(normals[i].data());
				assert(1.0 - n_i.norm() < FLT_EPSILON);
				ATA += n_i * n_i.transpose();
				assert(cube_bbox_eigen.contains(p_i));
				ATb += n_i * (n_i.transpose().dot(p_i));
			}

			Eigen::Vector3d cell_center = (cube_bbox_eigen.rbn + cube_bbox_eigen.ltf) / 2.0;
			Eigen::JacobiSVD<Eigen::Matrix3d> svd(ATA, Eigen::ComputeFullU | Eigen::ComputeFullV);
			svd.compute(ATA);
			Eigen::MatrixXd S = svd.singularValues().asDiagonal();
			// Invert the singular values diagonal matrix to calculate the pseudo inverse,
			// but cut off small singular values.
			for(short j = 0; j < 3; j++) {
				if(svd_cutoff_type == SVDCutoffType::RELATIVE && S(j, j) / S(0, 0) < svd_cutoff_factor) {
					// Relative cutoff like in Lindstrom 2000
					S(j, j) = 0.0;
				} else if(svd_cutoff_type == SVDCutoffType::ABSOLUTE && S(j, j) < svd_cutoff_factor) {
					// Absolute cutoff like in Ju et al. 2002
					S(j, j) = 0.0;
				} else {
					S(j, j) = 1.0 / S(j, j);
				}
			}
			Eigen::MatrixXd ATA_pinv = svd.matrixV() * S * svd.matrixU().transpose();
			Eigen::Vector3d x = cell_center + ATA_pinv * (ATb - ATA * cell_center);
			if(clip) {
				x = cube_bbox_eigen.clip(x, 0.2);
			}

			mesh.set_point(vh, OpenMesh::Vec3d(x.data()));
		}
#endif
	}
}