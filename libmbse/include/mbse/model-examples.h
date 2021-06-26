/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2021 University of Almeria                           |
  | Copyright (C) 2021 University of Salento                                |
  | See README for list of authors and papers                               |
  | Distributed under 3-clause BSD license                                  |
  |  See: <https://opensource.org/licenses/BSD-3-Clause>                    |
  +-------------------------------------------------------------------------+ */

#pragma once

#include <mbse/ModelDefinition.h>

namespace mbse
{
/** Builds a model of a NX x NY grid of four bar mechanisms.
 * \image docs/parameterized_fourbar.png
 */
ModelDefinition buildParameterizedMBS(
	const size_t nx, const size_t ny, const double NOISE_LEN = 0);

/** Builds a model of an N parts pendulum.
 * \todo Add picture!
 */
ModelDefinition buildLongStringMBS(
	const size_t N, double segmentLength = 0.5,
	double segmentMassPerMeter = 1.0);

/** Builds a model of an 4-bar mechanism.
 * \todo Add picture!
 * Degrees of freedom in q=[x1 y1 x2 y2]^T
 *
 */
ModelDefinition buildFourBarsMBS();

/** Builds a model of a slider crank mechanism.
 * \todo Add picture!
 * Degrees of freedom in q=[..]^T
 */
ModelDefinition buildSliderCrankMBS();

ModelDefinition buildFollowerMBS();

ModelDefinition buildTwoSliderBlocks();

}  // namespace mbse
