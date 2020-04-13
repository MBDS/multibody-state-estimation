/*+-------------------------------------------------------------------------+
  |    FactorGraph Control (sparsembs)  C++ library                         |
  |                                                                         |
  | Copyright (C) 2019-2020 University of Almeria                           |
  | See README for list of authors and papers                               |
  | Distributed under GNU General Public License version 3                  |
  |   See <http://www.gnu.org/licenses/>                                    |
  +-------------------------------------------------------------------------+ */

#pragma once

#include <mbse/CModelDefinition.h>

namespace mbse
{
/** Builds a model of a NX x NY grid of four bar mechanisms.
 * \image docs/parameterized_fourbar.png
 * \todo Add picture!
 */
void buildParameterizedMBS(
	const size_t nx, const size_t ny, CModelDefinition& model,
	const double NOISE_LEN = 0);

/** Builds a model of an N parts pendulum.
 * \todo Add picture!
 */
void buildLongStringMBS(const size_t N, CModelDefinition& model);

/** Builds a model of an 4-bar mechanism.
 * \todo Add picture!
 * Degrees of freedom in q=[x1 y1 x2 y2]^T
 *
 */
void buildFourBarsMBS(CModelDefinition& model);

/** Builds a model of a slider crank mechanism.
 * \todo Add picture!
 * Degrees of freedom in q=[..]^T
 */
void buildSliderCrankMBS(CModelDefinition& model);

void buildFollowerMBS(CModelDefinition& model);

void buildTwoSliderBlocks(CModelDefinition& model);

}  // namespace mbse
