#include "MishMesh/dual_contouring.h"

#include "vertexfit.hpp"
#include "grid.hpp"
#include <list>

namespace MishMesh {
	namespace Meshing {
		namespace DualContouring {

			struct MeshProperties {
				OpenMesh::VPropHandleT<std::vector<OpenMesh::Vec3d>> normals;
				OpenMesh::VPropHandleT<std::vector<OpenMesh::Vec3d>> intersection_points;
				OpenMesh::VPropHandleT<int> cube_idx;
			};

			/**
			 * An undirected edge
			 */
			struct Edge {
				int idx1, idx2;
				Edge(int idx1, int idx2) {
					this->idx1 = std::min(idx1, idx2);
					this->idx2 = std::max(idx1, idx2);
				}
				bool operator<(const Edge &other) const {
					return idx1 < other.idx1 || (idx1 == other.idx1 && idx2 < other.idx2);
				}
				bool operator==(const Edge &other) const {
					return idx1 != other.idx1 || idx2 != other.idx2;
				}
			};

			/**
			 * A directed halfedge storing everything needed to build a halfedge mesh.
			 */
			struct HalfEdge {
				int idx1, idx2;
				OpenMesh::Vec3d normal;
				OpenMesh::Vec3d intersection_point;
				HalfEdge *opposite;
				HalfEdge *next;
				HalfEdge *prev;
				MishMesh::PolyMesh::VertexHandle vh{};
				bool visited = false;
				bool boundary = false;
				bool quad_edge = false;
			};

			/**
			 * Test if an edge between points with these values is intersected by isosurface algorithms.
			 * This function checks for three conditions to resolve ambiguities:
			 * - If at least one point is infinite, there is no intersection.
			 * - If value1*value2 < 0, there is an intersection.
			 * - If one value is 0 and the other value is smaller than 0, there is an intersection.
			 * @param value1 The value of the first point.
			 * @param value2 The value of the second point.
			 * @returns True, if there is an intersection.
			 */
			inline bool has_intersection(double value1, double value2) {
				if(!std::isfinite(value1) || !std::isfinite(value2)) return false;
				if(value1 * value2 < 0) return true;
				if(value1 == 0 && value2 < 0) return true;
				if(value1 < 0 && value2 == 0) return true;
				return false;
			}

			/**
			 * Interpolate the point and normal along a grid edge from the coordinates and normals at the grid points.
			 * @param grid The grid.
			 * @param point_values The values at the grid points.
			 * @param point_normals The normals at the grid points..
			 * @param idx2 The index of the second grid point.
			 * @returns A interpolated point and the interpolated point normal.
			 */
			template<typename GridT>
			std::pair<OpenMesh::Vec3d, OpenMesh::Vec3d> interpolate_point_and_normal(const GridT &grid, const std::vector<double> &point_values, const std::vector<OpenMesh::Vec3d> &point_normals, int idx1, int idx2) {
				double value1 = point_values[idx1];
				double value2 = point_values[idx2];
				double alpha = -value2 / (value1 - value2);
				assert(alpha >= 0 && alpha <= 1);
				OpenMesh::Vec3d normal = (point_normals[idx1] * alpha + point_normals[idx2] * (1.0 - alpha)).normalized();
				OpenMesh::Vec3d point = grid.point(idx1) * alpha + grid.point(idx2) * (1.0 - alpha);
				return std::make_pair(point, normal);
			}

			/**
			 * Sample the distance function at the grid points.
<<<<<<< HEAD
			 * @param grid The grid.
			 * @param distance_function A function for calculating distances and normals.
>>>>>>> af85d32... Added a dual contouring implementation
			 */
			template<typename GridT>
			std::pair<std::vector<double>, std::vector<OpenMesh::Vec3d>> compute_point_values_and_normals(const GridT &grid, ::MishMesh::Meshing::DistanceFunction distance_function, void *distance_function_data) {
				const auto resolution = grid.resolution();
				std::vector<double> point_values(resolution[0] * resolution[1] * resolution[2]);
				std::vector<OpenMesh::Vec3d> point_normals(resolution[0] * resolution[1] * resolution[2]);
#pragma omp parallel for
				for(int i = 0; i < resolution[0]; i++) {
					for(int j = 0; j < resolution[1]; j++) {
						for(int k = 0; k < resolution[2]; k++) {
							OpenMesh::Vec3d p = grid.point(i, j, k);
							int index = grid.point_idx(i, j, k);
							std::tie(point_values[index], point_normals[index]) = distance_function(p, distance_function_data);
						}
					}
				}
				return make_pair(point_values, point_normals);
			}

			/**
			 * Create the halfedges belonging to one grid face, that are dual to the edges the face.
			 * @param[inout] halfedges The list with halfedges.
			 * @param[inout] gridedge_to_halfedges A map that assigns each grid edge the (four) halfedges that intersect
			 *               the faces adjacent to the grid edge.
			 * @param grid The grid.
			 * @param point_values The distance function values at the grid points.
			 * @param point_normals The distance function normals at the grid points.
			 * @param i The step in x direction.
			 * @param j The step in y direction.
			 * @param k The step in z direction.
			 * @param faceDirection The direction of the processed face. Each cube (i,j,k) is used to process its LEFT/TOP/FAR faces.
			 */
			template<typename GridT>
			void create_face_halfedges(std::list<HalfEdge> &halfedges, std::map<Edge, std::vector<HalfEdge *>> &gridedge_to_halfedges, const GridT &grid, const std::vector<double> &point_values, const std::vector<OpenMesh::Vec3d> &point_normals, int i, int j, int k, FaceDirection faceDirection) {
				// Get the indices of the corners of the face
				std::array<int, 4> idxs = grid.face_idxs(faceDirection, i, j, k);

				HalfEdge *previous_halfedge = nullptr;
				for(short l = 0; l < 4; l++) {
					// The edge intersects the isosurface
					if(has_intersection(point_values[idxs[l]], point_values[idxs[(l + 1) % 4]])) {
						OpenMesh::Vec3d point;
						OpenMesh::Vec3d normal;
						std::tie(point, normal) = interpolate_point_and_normal(grid, point_values, point_normals, idxs[l], idxs[(l + 1) % 4]);

						// Use an undirected primal cube edge to identify the quad around it.
						Edge edge = Edge(idxs[l], idxs[(l + 1) % 4]);

						// Create the halfedge between this cube and the left, top or far neighbor cube,
						// which has its index shifted by -1 in the corresponding dimension.
						int cube_idx1 = idxs[0];
						int cube_idx2 = grid.point_idx(cube_idx1,
						                               faceDirection == LEFT ? -1 : 0,
						                               faceDirection == TOP ? -1 : 0,
						                               faceDirection == FAR ? -1 : 0);

						// Flip the halfedge depending on the direction (positive to negative) of the primal cube edge.
						if(point_values[idxs[l]] > point_values[idxs[(l + 1) % 4]]) {
							std::swap(cube_idx1, cube_idx2);
						}

						halfedges.push_back({cube_idx1, cube_idx2,
						                     normal, point, previous_halfedge});

						auto halfedge_ptr = &halfedges.back();
						gridedge_to_halfedges[edge].push_back(halfedge_ptr);

						// A grid face is intersected by four dual halfedges, that belong to the grid face and one of its edges
						// and every second halfedge is connected to the previous one as its opposite edge.
						// This gives correct connectivity for the three possible cases:
						// 1) One different sign -> Create a wedge (e.g. crossing left and top edge)
						// 2) Two different signs next to each other -> Create plane (e.g. crossing left and right edge)
						// 3) Two different signs at opposite corners of the face -> Create two wedges
						//   (e.g. a wedge crossing left and top edges and a wedge crossing bottom and right edges)
						//
						//      |                         |
						//   +--|---+    +------+      +--|---+
						// __|__|   |  __|______|__  __|__|   |
						//   |      |    |      |      |    __|__
						//   |      |    |      |      |   |  |
						//   +------+    +------+      +---|--+
						//                                 |
						//     (1)        (2)             (3)
						if(previous_halfedge != nullptr) {
							assert(previous_halfedge->opposite == nullptr);
							previous_halfedge->opposite = halfedge_ptr;
							assert(halfedge_ptr->opposite == previous_halfedge);
							// Unset previous_halfedge, so the next intersection creates a new halfedge
							previous_halfedge = nullptr;
						} else {
							// Set previous_halfedge, so the next halfedge will be used as the opposite of this one.
							previous_halfedge = halfedge_ptr;
						}
					}
				}
			}

			/**
			 * Run create_face_halfedges for every face in the grid.
			 * @param[out] halfedges A list of halfedges.
			 * @param[out] gridedge_to_halfedge A map that contains for each grid edge a list of halfedges forming a face around the grid edge.
			 * @param grid The grid.
			 * @param point_values The distance values at the grid points.
			 * @param point_normals The normals at the grid points.
			 */
			template<typename GridT>
			void create_halfedges(std::list<HalfEdge> &halfedges, std::map<Edge, std::vector<HalfEdge *>> &gridedge_to_halfedges, const GridT &grid, const std::vector<double> &point_values, const std::vector<OpenMesh::Vec3d> &point_normals) {
				const auto resolution = grid.numCells();
				for(int i = 0; i < resolution[0]; i++) {
					for(int j = 0; j < resolution[1]; j++) {
						for(int k = 0; k < resolution[2]; k++) {
							for(short dim = 0; dim < 3; dim++) {
								FaceDirection faceDirection = dim == 0 ? FaceDirection::LEFT : (dim == 1 ? FaceDirection::TOP : FaceDirection::FAR);
								create_face_halfedges(halfedges, gridedge_to_halfedges, grid, point_values, point_normals, i, j, k, faceDirection);
							}
						}
					}
				}
			}

			/**
			 * Create quads in the halfedge structure using the information which halfedges belong to a dual face intersected by a grid edge.
			 * @param grid The grid.
			 * @param grid_edge_to_halfedges A map, that maps edges in the grid to the half edges that belong to the edge.
			 */
			template<typename GridT>
			void assemble_quads(const GridT &grid, const std::map<Edge, std::vector<HalfEdge *>> &grid_edge_to_halfedges) {
				for(auto &edge_quad_pair : grid_edge_to_halfedges) {
					auto &quad_halfedges = edge_quad_pair.second;
					if(quad_halfedges.size() != 4) {
						// Remove invalid quads, that are generated at the boundary
						// when a quad intersects cube faces that are outside of the bounding box
						for(auto &he : quad_halfedges) {
							he->quad_edge = false; // the edge does not belong to a quad
							if(he->opposite != nullptr) {
								// invalidate the opposite pointer of the opposite halfedge,
								// so that a boundary edge will be generated for the remaining half edge
								he->opposite->opposite = nullptr;
							}
						}
					} else {
						// assemble halfedges into a quad
						for(auto &he1 : quad_halfedges) {
							he1->quad_edge = true;
							for(auto &he2 : quad_halfedges) {
								if(he1->idx2 == he2->idx1) {
									he1->next = he2;
									he2->prev = he1;
								}
							}
						}
					}
				}
			}

			/**
			 * Create the boundary halfedges for the mesh, by creating a boundary edge for each halfedge that does not have an opposite edge
			 * and create the opposite, next and prev pointers.
			 * @param[inout] boundary_halfedges A list storing the created boundary halfedges.
			 * @param[inout] halfedges A list storing the inner halfedges.
			 * 				 The edges in the list will be changed with the correct adjacenct to the boundary half edges.
			 * @param grid The grid.
			 */
			template<typename GridT>
			void create_boundary_halfedges(std::list<HalfEdge> &boundary_halfedges, std::list<HalfEdge> &halfedges, const GridT &grid) {
				for(auto &he : halfedges) {
					if(he.opposite == nullptr) {
						boundary_halfedges.push_back({he.idx2, he.idx1, he.normal, he.intersection_point, &he});
						he.opposite = &boundary_halfedges.back();
						he.opposite->boundary = true;
					}
				}

				// Generate the next pointers for boundary halfedges by circulating around the start vertex
				// using the previous pointers of the inner halfedges.
				for(auto &he : boundary_halfedges) {
					HalfEdge *he2 = he.opposite;
					while(!he2->boundary) {
						he2 = he2->prev->opposite;
					}
					assert(he2->idx1 == he.idx2);
					he.next = he2;
					he2->prev = &he;
				}
			}

			/**
			 * For each vertex, collect the intersection points of adjacent edges and the normals at these points in vertex properties,
			 * so they can be used later on for calculating vertex positions.
<<<<<<< HEAD
			 * @param mesh The mesh.
			 * @param meshProperties A struct with the handles of the needed mesh properties.
			 * @param grid The grid.
			 * @param halfedges a list storing the halfedges.
>>>>>>> af85d32... Added a dual contouring implementation
			 */
			template<typename GridT>
			void collect_intersection_points_and_normals(MishMesh::PolyMesh &mesh, const MeshProperties &meshProperties, const GridT &grid, const std::list<HalfEdge> &halfedges) {
				for(auto &he : halfedges) {
					// Collect the interpolated normals and points from the edges to the vertex handle
					auto he2 = &he;
					// For each outgoing edge adjacent to the vertex that belongs to the halfedge he,
					// insert the intersection point and normal into the corresponding lists.
					do {
						he2 = he2->opposite->next;
						assert(he.vh.is_valid());
						mesh.property(meshProperties.normals, he.vh).push_back(he2->normal);
						mesh.property(meshProperties.intersection_points, he.vh).push_back(he2->intersection_point);
					} while(he2 != &he);
				}
			}

			/**
			 * Create vertices and faces from a halfedge data structure, by associating a new mesh vertex to each halfedge vertex, that does not
			 * have a vertex yet and creating faces for halfedge cycles.
			 * The vertices are initialized to lie on the center of the cuboid cells that contain them.
			 * @param mesh The mesh.
			 * @param meshProperties A struct with the handles of the needed mesh properties.
			 * @param grid The grid.
			 * @param halfedges a list storing the halfedges.
			 */
			template<typename GridT>
			void create_vertices_and_faces(MishMesh::PolyMesh &mesh, const MeshProperties &meshProperties, const GridT &grid, std::list<HalfEdge> &halfedges) {
				// We need to add all edges to the queue, because the input may have several disjunct connected components
				for(auto &he : halfedges) {
					if(!he.quad_edge) continue; // Ignore edges that are created at the boundary and do not belong to a full quad.
					if(he.visited) continue;    // Ignore halfedges that are already processed
					if(he.boundary) continue;   // Ignore boundary edges.
					auto quad_he = &he;         // Initialize the pointer to the current half edge.
					std::array<MishMesh::PolyMesh::VertexHandle, 4> vhs;
					for(short j = 0; j < 4; j++) {
						// If there is no vertex associated with the edge, try to find one on edges that share the same vertex.
						if(!quad_he->vh.is_valid()) {
							assert(quad_he->opposite != nullptr);
							assert(quad_he->opposite->next != nullptr);
							auto he2 = quad_he;
							do {
								// next edge in the 1-ring around the common vertex.
								assert(he2->opposite != nullptr);
								assert(he2->opposite->next != nullptr);
								he2 = he2->opposite->next;
								if(he2->vh.is_valid()) {
									quad_he->vh = he2->vh;
									break;
								}
							} while(he2 != quad_he);
						}
						// When no vertex was found, create one
						if(!quad_he->vh.is_valid()) {
							OpenMesh::Vec3i idxs_ltf = grid.grid_idxs(quad_he->idx1);
							OpenMesh::Vec3i idxs_rbn = idxs_ltf + OpenMesh::Vec3i{1, 1, 1};
							assert(idxs_rbn[0] < grid.resolution(0) && idxs_rbn[1] < grid.resolution(1) && idxs_rbn[2] < grid.resolution(2));
							auto p = grid.point(idxs_ltf) - grid.point(idxs_rbn);
							quad_he->vh = mesh.add_vertex(p);
							mesh.property(meshProperties.cube_idx, quad_he->vh) = quad_he->idx1;
						}
						vhs[j] = quad_he->vh;
						quad_he->visited = true;
						assert(quad_he->next != nullptr);
						quad_he = quad_he->next;
					}
					assert(quad_he == &he);
					mesh.add_face(vhs.data(), 4);
				}
				collect_intersection_points_and_normals(mesh, meshProperties, grid, halfedges);
			}

			/**
			 * Fit the created dual vertices to a position inside the grid cell.
			 * @param mesh The mesh.
			 * @param meshProperties A struct with the handles of the needed mesh properties.
			 * @param grid The grid.
			 * @param vertexFit The method for fitting the vertices.
			 */
			template<typename GridT>
			void fit_vertices(MishMesh::PolyMesh &mesh, const MeshProperties &meshProperties, const GridT &grid, const VertexFit vertexFit) {
				// Fit vertices to the position induced by the adjacent quads
				for(auto vh : mesh.vertices()) {
					auto idxs_ltf = grid.grid_idxs(mesh.property(meshProperties.cube_idx, vh));
					auto idxs_rbn = idxs_ltf + OpenMesh::Vec3i(1, 1, 1);
					MishMesh::BBox<OpenMesh::Vec3d, 3> vertex_bbox{grid.point(idxs_ltf), grid.point(idxs_rbn)};
					assert(vertex_bbox.is_valid());
					auto &points = mesh.property(meshProperties.intersection_points, vh);
					if(vertexFit == VertexFit::Average) {
						fit_vertex_average(mesh, vh, vertex_bbox, points);
					}
#ifdef HAS_EIGEN
					else if(vertexFit == VertexFit::SVD) {
						auto &normals = mesh.property(meshProperties.normals, vh);
						fit_vertex_svd(mesh, vh, vertex_bbox, points, normals, true, SVDCutoffType::ABSOLUTE, SVD_CUTOFF_FACTOR);
					}
#endif
				}
			}
		}

		/**
		 * Create a quadrangulation of the isosurface of a signed distance function using the dual contouring algorithm.
		 * @param distance_function A function that calculates the distance and the normal of a grid point.
		 * @param distance_function_data An optional pointer to data used by the distance_function.
		 * @param mesh_bounding_box The bounding box of the mesh.
		 * @param resolution The number of grid cells in each dimension.
		 * @param vertexFit The method for fitting the vertices.
		 */
		MishMesh::PolyMesh dual_contouring(DistanceFunction distance_function, void *distance_function_data, MishMesh::BBox<OpenMesh::Vec3d, 3> mesh_bounding_box, OpenMesh::Vec3i resolution, VertexFit vertexFit) {
			auto grid = RegularGrid(mesh_bounding_box, resolution + OpenMesh::Vec3i{1, 1, 1});
			MishMesh::PolyMesh mesh;

			DualContouring::MeshProperties meshProperties;
			mesh.add_property(meshProperties.cube_idx);
			mesh.add_property(meshProperties.intersection_points);
			mesh.add_property(meshProperties.normals);

			// Calculate the values and normals at the grid points
			std::vector<double> point_values;
			std::vector<OpenMesh::Vec3d> point_normals;
			std::tie(point_values, point_normals) = DualContouring::compute_point_values_and_normals(grid, distance_function, distance_function_data);

			// Create dual halfedges around edges that intersect the zero-contour
			std::list<DualContouring::HalfEdge> halfedges;
			std::map<DualContouring::Edge, std::vector<DualContouring::HalfEdge *>> gridedge_to_halfedges;
			DualContouring::create_halfedges(halfedges, gridedge_to_halfedges, grid, point_values, point_normals);

			// Assemble the halfedges to quads
			DualContouring::assemble_quads(grid, gridedge_to_halfedges);

			// Remove halfedges that do not belong to a quad
			halfedges.remove_if([](const DualContouring::HalfEdge &he) { return !(he.quad_edge); });

			// Create mesh boundary halfedges
			std::list<DualContouring::HalfEdge> boundary_halfedges;
			DualContouring::create_boundary_halfedges(boundary_halfedges, halfedges, grid);

			for(auto &he : halfedges) {
				assert(he.opposite != nullptr);
			}

			// Create the mesh vertices and faces using the halfedge structure
			DualContouring::create_vertices_and_faces(mesh, meshProperties, grid, halfedges);

			// Fit the vertices to their positions, either averaged or fit to a plane
			DualContouring::fit_vertices(mesh, meshProperties, grid, vertexFit);

			return mesh;
		}
	}
}
