/*
 * Copyright (C) 2024-2025, Alexander Schier
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#ifndef ACC_TRAITS_HEADER
#define ACC_TRAITS_HEADER

#include <type_traits>
#include "defines.h"

ACC_NAMESPACE_BEGIN

// Types with square_norm like math/vector.h
template<typename T, typename = void>
struct has_square_norm : std::false_type {};
template<typename T>
struct has_square_norm<T, std::void_t<decltype(std::declval<T>().square_norm())>> : std::true_type {};

// Types with squaredNorm like Eigen::Vector3d
template<typename T, typename = void>
struct has_squaredNorm : std::false_type {};
template<typename T>
struct has_squaredNorm<T, std::void_t<decltype(std::declval<T>().squaredNorm())>> : std::true_type {};

// Types with sqrnorm like OpenMesh::Vec3d
template<typename T, typename = void>
struct has_sqrnorm : std::false_type {};
template<typename T>
struct has_sqrnorm<T, std::void_t<decltype(std::declval<T>().sqrnorm())>> : std::true_type {};

// Types that are unsupported
template<typename>
struct dependent_false : std::false_type {};


template<typename T>
double squaredNorm(const T& vec) {
	if constexpr (has_square_norm<T>::value) {
		return vec.square_norm();
	} else if constexpr (has_squaredNorm<T>::value) {
		return vec.squaredNorm();
	} else if constexpr (has_sqrnorm<T>::value) {
		return vec.sqrnorm();
	} else {
		static_assert(dependent_false<T>::value, "Norm function not implemented for this type");
		return 0.0;
	}
}

ACC_NAMESPACE_END

#endif /* ACC_TRAITS_HEADER */
