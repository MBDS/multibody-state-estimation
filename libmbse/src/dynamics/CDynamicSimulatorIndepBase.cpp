/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2020 University of Almeria                           |
  | Copyright (C) 2020 University of Salento                                |
  | See README for list of authors and papers                               |
  | Distributed under 3-clause BSD license                                  |
  |  See: <https://opensource.org/licenses/BSD-3-Clause>                    |
  +-------------------------------------------------------------------------+ */

#include <mbse/CModelDefinition.h>
#include <mbse/CAssembledRigidModel.h>
#include <mbse/dynamics/dynamic-simulators.h>

using namespace mbse;
using namespace Eigen;
using namespace mrpt::math;
using namespace mrpt;
using namespace std;

// Ctor
CDynamicSimulatorIndepBase::CDynamicSimulatorIndepBase(
	std::shared_ptr<CAssembledRigidModel> arm_ptr)
	: CDynamicSimulatorBase(arm_ptr)
{
}

// Virtual destructor: requird for virtual bases
CDynamicSimulatorIndepBase::~CDynamicSimulatorIndepBase() {}

void CDynamicSimulatorIndepBase::solve_ddotz(double t, VectorXd& ddot_z)
{
	ASSERT_(init_);
	this->internal_solve_ddotz(t, ddot_z);
}

// Run simulation:
double CDynamicSimulatorIndepBase::run(const double t_ini, const double t_end)
{
	ASSERT_(t_end >= t_ini);
	ASSERT_(init_);

	if (t_ini == t_end) return t_end;  // Nothing to do.

	// Fill constant data
	TSimulationState sim_state(arm_);

	const double t_step = std::min(t_end - t_ini, params.time_step);
	const double t_step2 = t_step * 0.5;
	const double t_step6 = t_step / 6.0;

	double t;  // Declared here so we know the final "time":
	for (t = t_ini; t < t_end; t += t_step)
	{
		// Log sensor points:
		// ------------------------------
		for (std::list<TSensorData>::iterator it = sensors_.begin();
			 it != sensors_.end(); ++it)
		{
			TSensorData& sd = *it;
			sd.log.push_back(timestamped_point_t(
				t, TPointState(
					   mrpt::math::TPoint2D(*sd.pos[0], *sd.pos[1]),
					   mrpt::math::TPoint2D(*sd.vel[0], *sd.vel[1]))));
		}

		timelog().enter("mbs.run_complete_timestep");

		// Integrate:
		// ------------------------------
		switch (params.ode_solver)
		{
			case ODE_Euler:
			{
				this->internal_solve_ddotz(t, ddotz1);
				arm_->q_ += t_step * arm_->dotq_;
				// arm_->dotq_ += t_step * ddotz1;
				this->dq_plus_dz(arm_->dotq_, t_step * ddotz1, arm_->dotq_);

				this->correct_dependent_q_dq();
			}
			break;

			case ODE_RK4:
			{
				q0 = arm_->q_;  // Make backup copy of state (velocities will
								// be in "v1")

				// k1 = f(t,y);
				// cur_time = t;
				v1 = arm_->dotq_;
				// No change needed: arm_->q_ = q0;
				can_choose_indep_coords_ = true;
				this->internal_solve_ddotz(t, ddotz1);
				can_choose_indep_coords_ = false;  // don't change indep. coords

				// k2 = f(t+At/2,y+At/2*k1)
				// cur_time = t + t_step2;
				this->dq_plus_dz(
					v1, t_step2 * ddotz1,
					arm_->dotq_);  // \dot{q}= \dot{q}_0 + At/2 * \ddot{q}_1
				arm_->q_ = q0 + t_step2 * v1;
				this->correct_dependent_q_dq();

				v2 = arm_->dotq_;
				this->internal_solve_ddotz(t + t_step2, ddotz2);

				// k3 = f(t+At/2,y+At/2*k2)
				// cur_time = t + t_step2;
				// arm_->dotq_ = v1 + t_step2*ddotq2;  // \dot{q}= \dot{q}_0 +
				// At/2 * \ddot{q}_2
				this->dq_plus_dz(v1, t_step2 * ddotz2, arm_->dotq_);

				arm_->q_ = q0 + t_step2 * v2;
				this->correct_dependent_q_dq();

				v3 = arm_->dotq_;
				this->internal_solve_ddotz(t + t_step2, ddotz3);

				// k4 = f(t+At  ,y+At*k3)
				// cur_time = t + t_step;
				// arm_->dotq_ = v1 + t_step*ddotq3;
				this->dq_plus_dz(v1, t_step * ddotz3, arm_->dotq_);
				arm_->q_ = q0 + t_step * v3;
				this->correct_dependent_q_dq();

				v4 = arm_->dotq_;
				this->internal_solve_ddotz(t + t_step, ddotz4);

				// Runge-Kutta 4th order formula:
				arm_->q_ = q0 + t_step6 * (v1 + 2 * v2 + 2 * v3 + v4);
				this->dq_plus_dz(
					v1, t_step6 * (ddotz1 + 2 * ddotz2 + 2 * ddotz3 + ddotz4),
					arm_->dotq_);
				this->correct_dependent_q_dq();
			}
			break;

			// Implicit trapezoidal integration rule:
			// -------------------------------------------
#if 0
		case ODE_Trapezoidal:
			{
				const double t_step_sq = t_step * t_step;

				const size_t MAX_ITERS = 10;
				const double QDIFF_MAX = 1e-10;
				double qdiff=10*QDIFF_MAX;

				// Keep the initial state:
				const Eigen::VectorXd q0  = arm_->q_;
				const Eigen::VectorXd dq0 = arm_->dotq_;

				// First attempt:
				Eigen::VectorXd ddz0;
				this->can_choose_indep_coords_=true;
				this->internal_solve_ddotz(ddz0, true);
				this->can_choose_indep_coords_=false;

				Eigen::VectorXd  q_new =  q0 + t_step*dq0 + 0.5*t_step_sq*ddq0;
				Eigen::VectorXd dq_new;
				this->dq_plus_dz(dq0, t_step*ddq0, dq_new );

				// Solve at the new predicted state "t=k+1":
				arm_->q_    =  q_new;
				arm_->dotq_ = dq_new;
				this->correct_dependent_q_dq();
				Eigen::VectorXd q_old = q_new;

				Eigen::VectorXd ddq_mid;
				size_t iter;
				for (iter=0; iter<MAX_ITERS && qdiff>QDIFF_MAX ; iter++ )
				{
					// Store previous state for comparing the progress of the iterative method:
					q_old = q_new;

					// Solve at the new predicted state "t=k+1":
					this->internal_solve_ddotz(ddotz1);
					//this->solve_ddotq(ddotq1);

					// integrator (trapezoidal rule)
					// -------------------------------
					ddq_mid = (ddotz1+ddz0)*0.5;
					q_new =  q0 + t_step*dq0 + 0.5*t_step_sq*ddq_mid;
					dq_new = dq0 + t_step*ddq_mid;

					// check progress:
					qdiff = (q_old-q_new).norm();

					// Solve at the new predicted state "t=k+1":
					arm_->q_    =  q_new;
					arm_->dotq_ = dq_new;
				}
				ASSERTMSG_(iter<MAX_ITERS,"Trapezoidal convergence failed!")

				timelog().registerUserMeasure("trapezoidal.iters",iter);
			}
			break;
#endif

			default:
				THROW_EXCEPTION("Unknown value for params.ode_solver");
		};

		// Save last ddotq:
		CAssembledRigidModel::TComputeDependentResults cdr;
		cdr.ddotq = &arm_->ddotq_;
		arm_->computeDependentPosVelAcc(
			independent_coordinate_indices(), false /*update q*/,
			false /*update dq*/, {}, cdr, &ddotz1);

		timelog().leave("mbs.run_complete_timestep");

		// User-callback:
		// ------------------------------
		sim_state.t = t;
		if (params.user_callback) params.user_callback(sim_state);
	}

	return t;
}

/** Wrapper for ddotq computation, from ddotz */
void CDynamicSimulatorIndepBase::internal_solve_ddotq(
	double t, VectorXd& ddot_q, VectorXd* lagrangre)
{
	ASSERTMSG_(
		lagrangre == nullptr,
		"This solver uses independent coordinates, so it can't determine the "
		"lagrange multipliers as requested!");

	THROW_EXCEPTION("TO DO! Better use internal_solve_ddotz() instead.");
}
