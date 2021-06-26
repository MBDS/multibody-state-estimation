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

#include <gtsam/base/Vector.h>
#include <gtsam/base/VectorSpace.h>

namespace mbse
{
/** Type for system internal states q_{k}, dq_{k}, ddq_{k} */
using state_t = gtsam::Vector;

}  // namespace mbse
