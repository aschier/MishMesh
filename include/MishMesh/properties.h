#pragma once

#include <string>
#include <OpenMesh/Core/Utils/Property.hh>

namespace MishMesh {
	/**
	 * Get a property handle with a given name from a mesh.
	 * @param mesh The mesh.
	 * @param property_name The name of the property.
	 * @tparam PropertyHandleT The type of the property.
	 * @tparam MeshT The mesh type. Can usually be deducted automatically.
	 */
	template<typename PropertyHandleT, typename MeshT>
	inline PropertyHandleT get_property_handle(const MeshT &mesh, const std::string &property_name) {
		PropertyHandleT property;
		mesh.get_property_handle(property, property_name);
		return property;
	}

	/**
	 * Add a property with a given name to a mesh and returns the created property.
	 *
	 * @param mesh The mesh.
	 * @param property_name The name of the property.
	 * @tparam PropertyT The type of the property.
	 * @tparam MeshT The mesh type. Can usually be deduced automatically.
	 * @return PropertyT& The created property.
	 */
	template<typename PropertyT, typename MeshT>
	inline PropertyT add_property(MeshT &mesh, const std::string &property_name) {
		PropertyT property;
		mesh.add_property(property, property_name);
		return property;
	}

	/**
	 * Copy a face property from a mesh to a submesh of the mesh.
	 * @param base_mesh The base mesh
	 * @param[inout] sub_mesh A submesh of the base_mesh with matching face ids.
	 * @param submesh_to_basemesh_property A face property that stores the base_mesh indices on vertices of the submesh.
	 * @param basemesh_property The property to copy from the base mesh.
	 * @param submesh_property The property on the submesh.
	 * @tparam The mesh type
	 * @tparam The type of the property to copy
	 */
	template<typename MeshT, typename T>
	inline void copy_property(const MeshT &base_mesh, MeshT &sub_mesh, const OpenMesh::FPropHandleT<T> submesh_to_basemesh_property,
	                          const OpenMesh::FPropHandleT<T> base_mesh_property, const OpenMesh::FPropHandleT<int> submesh_property) {
		for(const auto &sub_mesh_fh : sub_mesh.faces()) {
			int base_mesh_idx = sub_mesh.property(submesh_to_basemesh_property, sub_mesh_fh);
			const auto &base_mesh_fh = base_mesh.face_handle(base_mesh_idx);
			sub_mesh.property(submesh_property, sub_mesh_fh) = base_mesh.property(base_mesh_property, base_mesh_fh);
		}
	}

	/**
	 * Copy a vertex property from a mesh to a submesh of the mesh.
	 * @param base_mesh The base mesh
	 * @param[inout] sub_mesh A submesh of the base_mesh with matching vertex ids.
	 * @param submesh_to_basemesh_property A vertex property that stores the base_mesh indices on vertices of the submesh.
	 * @param basemesh_property The property to copy from the base mesh.
	 * @param submesh_property The property on the submesh.
	 * @tparam The mesh type
	 * @tparam The type of the property to copy
	 */
	template<typename MeshT, typename T>
	inline void copy_property(const MeshT &base_mesh, MeshT &sub_mesh, const OpenMesh::VPropHandleT<T> submesh_to_basemesh_property,
	                          const OpenMesh::VPropHandleT<T> base_mesh_property, const OpenMesh::VPropHandleT<int> submesh_property) {
		for(const auto &sub_mesh_vh : sub_mesh.vertices()) {
			int base_mesh_idx = sub_mesh.property(submesh_to_basemesh_property, sub_mesh_vh);
			const auto &base_mesh_vh = base_mesh.vertex_handle(base_mesh_idx);
			sub_mesh.property(submesh_property, sub_mesh_vh) = base_mesh.property(base_mesh_property, base_mesh_vh);
		}
	}

	/**
	 * Requests the status property for all types of elements in the mesh.
	 *
	 * @tparam MeshT The mesh type.
	 * @param mesh The mesh.
	 */
	template<typename MeshT>
	inline void request_all_status(MeshT &mesh) {
		mesh.request_vertex_status();
		mesh.request_halfedge_status();
		mesh.request_edge_status();
		mesh.request_face_status();
	}

	/**
	 * Releases the status property for all types of elements in the mesh.
	 *
	 * @tparam MeshT The mesh type.
	 * @param mesh The mesh.
	 */
	template<typename MeshT>
	inline void release_all_status(MeshT &mesh) {
		mesh.release_vertex_status();
		mesh.release_halfedge_status();
		mesh.release_edge_status();
		mesh.release_face_status();
	}
}
