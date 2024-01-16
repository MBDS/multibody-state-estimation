/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2024 University of Almeria                           |
  | Copyright (C) 2021 University of Salento                                |
  | See README for list of authors and papers                               |
  | Distributed under 3-clause BSD license                                  |
  |  See: <https://opensource.org/licenses/BSD-3-Clause>                    |
  +-------------------------------------------------------------------------+ */

#include <gtest/gtest.h>

#include <mbse/model-examples.h>
#include <mbse/dynamics/dynamic-simulators.h>
#include <mbse/AssembledRigidModel.h>
#include <gtsam/inference/Symbol.h>
#include <gtsam/nonlinear/factorTesting.h>
#include <mbse/factors/FactorConstraintsVelIndep.h>

using namespace std;
using namespace mbse;

MRPT_TODO(
	"Enhance factor tests to also run with other mechanism types with more "
	"constraint types.");

TEST(Jacobians, FactorVelConstraintsIndep)
{
	MRPT_START

	// for use in EXPECT_CORRECT_FACTOR_JACOBIANS
	const auto name_ = "FactorVelConstraintsIndep";
#define EXPECT ASSERT_

	const auto Q = gtsam::SymbolGenerator('q');
	const auto dotQ = gtsam::SymbolGenerator('v');
	// const auto ddotQ = gtsam::SymbolGenerator('a');
	// const auto Z = gtsam::SymbolGenerator('z');
	const auto dotZ = gtsam::SymbolGenerator('x');
	// const auto ddotZ = gtsam::SymbolGenerator('c');
	using namespace mbse;

	// Create the multibody object:
	ModelDefinition model = mbse::buildFourBarsMBS();

	// Add an extra relative coordinate:
	model.rDOFs_.emplace_back(mbse::RelativeAngleAbsoluteDOF(0, 1));

	auto aMBS = model.assembleRigidMBS();
	aMBS->setGravityVector(0, -9.81, 0);

	CDynamicSimulator_Indep_dense dynSimul(aMBS);

	// Enforce the use of the theta angle as independent coordinate:
	dynSimul.independent_coordinate_indices({4});

	// Must be called before solve_ddotq():
	dynSimul.prepare();

	for (int ti = 0; ti < 10; ti++)
	{
		const double dt = 1.0;
		double t = ti * dt;
		dynSimul.run(t, t + dt);

		std::cout << "Evaluating test for t=" << t << "\n";
		std::cout << "q  =" << aMBS->q_.transpose() << "\n";
		std::cout << "dq =" << aMBS->dotq_.transpose() << "\n";
		std::cout << "ddq =" << aMBS->ddotq_.transpose() << "\n";

		// Add factors:
		// Create factor noises:
		// const auto n = aMBS->q_.size();
		const auto m = aMBS->Phi_q_.getNumRows();
		const std::vector<size_t> indCoordsInd = {4};
		// indCoordsInd = dynSimul.independent_coordinate_indices();
		const auto d = indCoordsInd.size();
		EXPECT_GT(d, 0);

		auto noise = gtsam::noiseModel::Isotropic::Sigma(m + d, 0.1);

		// Create a dummy factor:
		const auto factor = FactorConstraintsVelIndep(
			aMBS, indCoordsInd, noise, Q(1), dotQ(1), dotZ(1));

		// Convert plain Eigen vectors into state_t classes (used as Values
		// in GTSAM factor graphs):
		const state_t q = state_t(aMBS->q_);
		const state_t dotq = state_t(aMBS->dotq_);
		const state_t dotz = state_t(mbse::subset(dotq, indCoordsInd));

		gtsam::Values values;
		values.insert(Q(1), q);
		values.insert(dotQ(1), dotq);
		values.insert(dotZ(1), dotz);

		EXPECT_CORRECT_FACTOR_JACOBIANS(factor, values, 1e-9, 1e-3);
	}
	MRPT_END
}
