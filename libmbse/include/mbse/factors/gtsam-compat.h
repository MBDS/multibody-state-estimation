/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2024 University of Almeria                           |
  | Copyright (C) 2021 University of Salento                                |
  | See README for list of authors and papers                               |
  | Distributed under 3-clause BSD license                                  |
  |  See: <https://opensource.org/licenses/BSD-3-Clause>                    |
  +-------------------------------------------------------------------------+ */

#pragma once

#include <gtsam/config.h>

/** True if the active GTSAM version is >= major.minor.patch. */
#define GTSAM_VERSION_AT_LEAST(major, minor, patch)                       \
	((GTSAM_VERSION_MAJOR > (major)) ||                                   \
	 (GTSAM_VERSION_MAJOR == (major) && GTSAM_VERSION_MINOR > (minor)) || \
	 (GTSAM_VERSION_MAJOR == (major) && GTSAM_VERSION_MINOR == (minor) && \
	  GTSAM_VERSION_PATCH >= (patch)))

/** GTSAM < 4.3 declares the optional Jacobian outputs of NoiseModelFactorN
 * (aliased here as NoiseModelFactor1..6) as boost::optional<Matrix&>; GTSAM
 * >= 4.3 replaced that with a plain Matrix* (OptionalMatrixType/OptionalNone).
 * MBSE_OptionalMatrixType follows whichever the active GTSAM expects, so a
 * single evaluateError() override compiles against both without changes.
 * Since boost::optional<Matrix&> and Matrix* share the same `if(H)` / `*H`
 * usage, factor bodies never need to branch on this macro themselves.
 */
#define MBSE_GTSAM_USES_BOOST_OPTIONAL (!GTSAM_VERSION_AT_LEAST(4, 3, 0))

#if MBSE_GTSAM_USES_BOOST_OPTIONAL
#include <boost/optional.hpp>
#include <boost/pointer_cast.hpp>
#define MBSE_OptionalMatrixType boost::optional<gtsam::Matrix&>
#ifndef OptionalNone
#define OptionalNone boost::none
#endif
#else
// The NoiseModelFactor1..6 aliases and OptionalMatrixType/OptionalNone are
// defined either directly in NonlinearFactor.h, or (some GTSAM >=4.3
// versions) in the separate header below, which not every such version
// ships -- include it only if present.
#include <gtsam/nonlinear/NonlinearFactor.h>
#if defined(__has_include)
#if __has_include(<gtsam/nonlinear/NoiseModelFactorN.h>)
#include <gtsam/nonlinear/NoiseModelFactorN.h>
#endif
#else
#include <gtsam/nonlinear/NoiseModelFactorN.h>
#endif
#define MBSE_OptionalMatrixType gtsam::OptionalMatrixType
#endif

/** GTSAM < 4.3 uses boost::shared_ptr for gtsam::NonlinearFactor::shared_ptr
 * (and thus for This::clone()'s return value); GTSAM >= 4.3 uses
 * std::shared_ptr. MBSE_static_pointer_cast<T>(x) follows whichever the
 * active GTSAM expects. */
#if MBSE_GTSAM_USES_BOOST_OPTIONAL
#define MBSE_static_pointer_cast boost::static_pointer_cast
#else
#define MBSE_static_pointer_cast std::static_pointer_cast
#endif

/** Wraps a `gtsam::Matrix` lvalue as whatever evaluateError() expects for a
 * Jacobian output argument at a call site: a plain reference under GTSAM
 * <4.3 (implicitly convertible to boost::optional<Matrix&>), or its address
 * under GTSAM >=4.3 (Matrix*). */
#if MBSE_GTSAM_USES_BOOST_OPTIONAL
#define MBSE_JACARG(x) (x)
#else
#define MBSE_JACARG(x) (&(x))
#endif
