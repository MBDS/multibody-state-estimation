
#pragma once

#include <sparsembs/factor-common.h>
#include <sparsembs/dynamic-simulators.h>
#include <gtsam/nonlinear/NonlinearFactor.h>

namespace sparsembs
{
/** Factor for contraints
 * Here the constraints equation Phi(q)=0 is implmented
 */

// Create derived class "ConstraintsFactor" from superclass "NoiseModelFactor4"
class ConstraintsFactor
    : public gtsam::NoiseModelFactor4<state_t, state_t, state_t, state_t>
{
};
};  // namespace sparsembs
