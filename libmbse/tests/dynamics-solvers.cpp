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

#include <mbse/mbse.h>
#include <mbse/model-examples.h>

template <class DYNAMIC_SOLVER_T>
void testerPendulumDynamics(bool addRelativeAngle = false)
{
	const double L = 0.5;  // [m]
	const double massPerL = 1.0;  // [kg/m]

	mbse::timelog().enable(false);	// avois clutter in cout

	mbse::ModelDefinition model = mbse::buildLongStringMBS(1, L, massPerL);

	// optional relative DOFs:
	if (addRelativeAngle)
	{
		// Add angle between points #0 and #1:
		model.rDOFs_.emplace_back(mbse::RelativeAngleAbsoluteDOF(0, 1));
	}

	std::shared_ptr<mbse::AssembledRigidModel> aMBS = model.assembleRigidMBS();

	aMBS->setGravityVector(0, -9.81, 0);

	DYNAMIC_SOLVER_T dynSimul(aMBS);

	dynSimul.prepare();

	Eigen::VectorXd ddotq0;
	dynSimul.solve_ddotq(0.0 /*current time*/, ddotq0);

	const Eigen::VectorXd ddotq_real =
		addRelativeAngle ?	//
			Eigen::VectorXd(
				(Eigen::Vector3d() << 0, -14.7041, -29.43).finished())
						 :	//
			Eigen::VectorXd((Eigen::Vector2d() << 0, -14.7041).finished());

	EXPECT_NEAR(
		(ddotq0 - ddotq_real).array().abs().maxCoeff() /
			ddotq_real.array().abs().maxCoeff(),
		0, 1e-3)
		<< "ddotq0     : " << ddotq0.transpose() << "\n"
		<< "ddotq_real : " << ddotq_real.transpose() << "\n";
}

TEST(PendulumDynamics, CDynamicSimulator_Lagrange_LU_dense)
{
	testerPendulumDynamics<mbse::CDynamicSimulator_Lagrange_LU_dense>();
}
TEST(PendulumDynamics, CDynamicSimulator_Lagrange_UMFPACK)
{
	testerPendulumDynamics<mbse::CDynamicSimulator_Lagrange_UMFPACK>();
}
TEST(PendulumDynamics, CDynamicSimulator_Lagrange_KLU)
{
	testerPendulumDynamics<mbse::CDynamicSimulator_Lagrange_KLU>();
}
TEST(PendulumDynamics, CDynamicSimulator_Lagrange_CHOLMOD)
{
	testerPendulumDynamics<mbse::CDynamicSimulator_Lagrange_CHOLMOD>();
}
TEST(PendulumDynamics, CDynamicSimulator_AugmentedLagrangian_KLU)
{
	testerPendulumDynamics<mbse::CDynamicSimulator_AugmentedLagrangian_KLU>();
}
TEST(PendulumDynamics, CDynamicSimulator_AugmentedLagrangian_Dense)
{
	testerPendulumDynamics<mbse::CDynamicSimulator_AugmentedLagrangian_Dense>();
}
TEST(PendulumDynamics, CDynamicSimulator_ALi3_Dense)
{
	testerPendulumDynamics<mbse::CDynamicSimulator_ALi3_Dense>();
}
TEST(PendulumDynamics, CDynamicSimulator_R_matrix_dense)
{
	testerPendulumDynamics<mbse::CDynamicSimulator_R_matrix_dense>();
}

// ---------
TEST(PendulumDynamicsWithRelCoord, CDynamicSimulator_Lagrange_LU_dense)
{
	testerPendulumDynamics<mbse::CDynamicSimulator_Lagrange_LU_dense>(true);
}
TEST(PendulumDynamicsWithRelCoord, CDynamicSimulator_R_matrix_dense)
{
	testerPendulumDynamics<mbse::CDynamicSimulator_R_matrix_dense>(true);
}

// -------------
// test that acceleration from dynamics doesn't suffer from any discontinuity
// at the change of formulation for the absolute angle constraint
// (at sin(th)=0.707)
TEST(PendulumDynamicsWithRelCoord, NoDiscontinuityAcceleration)
{
	const double L = 0.5;  // [m]
	const double massPerL = 1.0;  // [kg/m]

	mbse::timelog().enable(false);	// avois clutter in cout

	mbse::ModelDefinition model = mbse::buildLongStringMBS(1, L, massPerL);

	// Add angle between points #0 and #1:
	model.rDOFs_.emplace_back(mbse::RelativeAngleAbsoluteDOF(0, 1));

	std::shared_ptr<mbse::AssembledRigidModel> aMBS = model.assembleRigidMBS();

	aMBS->setGravityVector(0, -9.81, 0);

	const std::vector<double> thetas = {
		asin(0.707) - 0.0015, asin(0.707) - 0.0005, asin(0.707) + 0.0005,
		asin(0.707) + 0.0015};

	std::vector<Eigen::VectorXd> ddqs;

	for (const double theta : thetas)
	{
		aMBS->q_(2) = theta;
		aMBS->ddotq_(2) = 1.0;
		mbse::AssembledRigidModel::ComputeDependentParams cdp;
		mbse::AssembledRigidModel::ComputeDependentResults cdr;
		Eigen::VectorXd ddotz, ddotq;
		ddotz.setConstant(1, 1.0);
		cdr.ddotq = &ddotq;
		aMBS->computeDependentPosVelAcc({2}, true, true, cdp, cdr, &ddotz);

		mbse::CDynamicSimulator_R_matrix_dense dynSimul(aMBS);
		dynSimul.prepare();

		Eigen::VectorXd ddotq0;
		dynSimul.solve_ddotq(0.0 /*current time*/, ddotq0);

#if 0
		std::cout << "q            : " << aMBS->q_.transpose() << std::endl;
		std::cout << "ddq(forced)  : " << cdr.ddotq->transpose() << std::endl;
		std::cout << "ddq(dynamics): " << ddotq0.transpose() << std::endl;
#endif
		ddqs.push_back(ddotq0);
	}

	const Eigen::VectorXd d0 = ddqs[1] - ddqs[0];
	const Eigen::VectorXd d1 = ddqs[2] - ddqs[1];
	const Eigen::VectorXd d2 = ddqs[3] - ddqs[2];

	EXPECT_NEAR(d0(2), d1(2), 1e-3);
	EXPECT_NEAR(d1(2), d2(2), 1e-3);
}
