/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2021 University of Almeria                           |
  | Copyright (C) 2021 University of Salento                                |
  | See README for list of authors and papers                               |
  | Distributed under 3-clause BSD license                                  |
  |  See: <https://opensource.org/licenses/BSD-3-Clause>                    |
  +-------------------------------------------------------------------------+ */

#include <mbse/MultiBodyParticleFilter.h>

#include <mrpt/math/distributions.h>

using namespace mbse;
using namespace Eigen;
using namespace mrpt::math;
using namespace mrpt;
using namespace std;

// Ctor:
MultiBodyParticleFilter::MultiBodyParticleFilter(
	const size_t M, const ModelDefinition& mbs)
{
	// 1) Proccess model:
	TSymbolicAssembledModel sym_model(mbs);
	mbs.assembleRigidMBS(sym_model);

	// 2) Create particles:
	m_particles.resize(M);

	for (auto& p : m_particles)
	{
		p.log_w = 0;
		p.d.reset(new particle_t(sym_model));
	}

	// Randomize?
	random_generator.randomize();
}

// Dtor:
MultiBodyParticleFilter::~MultiBodyParticleFilter() {}

void MultiBodyParticleFilter::run_PF_step(
	const double t_ini, const double t_end, const double max_t_step,
	const std::vector<CVirtualSensor::Ptr>& sensor_descriptions,
	const std::vector<double>& sensor_readings, TOutputInfo& out_info)
{
	mrpt::system::CTimeLoggerEntry tle(timelog(), "run_PF_step");

	ASSERT_(sensor_descriptions.size() == sensor_readings.size());

	// 1) Executes probabilistic transition model:
	// -----------------------------------------------------
	timelog().enter("PF.1.forward_model");

	ASSERT_ABOVE_(t_end, t_ini);
	const double t_increment = t_end - t_ini;
	const size_t nTimeSteps = ceil(t_increment / max_t_step);

	const double t_step = t_increment / nTimeSteps;
	const double t_step2 = t_step * 0.5;
	const double t_step6 = t_step / 6.0;

	Eigen::VectorXd q0;	 // Backup of state.
	Eigen::VectorXd k1, k2, k3, k4;
	Eigen::VectorXd v1, v2, v3, v4;	 // \dot{q}
	Eigen::VectorXd ddotz1, ddotz2, ddotz3, ddotz4;	 // \ddot{q}
	Eigen::VectorXd q_incr, dotz_incr;
	Eigen::VectorXd dotz_noise;

	const size_t nParts = m_particles.size();
	size_t i;
	double t = t_ini;
	for (size_t nTim = 0; nTim < nTimeSteps; nTim++, t += t_step)
	{
#if 0 && MBS_HAVE_OPENMP
#pragma omp                                                                                                               \
	parallel for private(q0) private(k1) private(k2) private(k3) private(k4) private(v1) private(v2) private(v3) private( \
		v4) private(ddotq1) private(ddotq2) private(ddotq3) private(ddotq4) private(q_incr) private(dotq_incr) private(dotq_noise)
#endif
		for (i = 0; i < nParts; i++)
		{
			auto part = m_particles[i].d;

			// ODE_RK4:
			// --------------------------------
			{
				q0 = part->num_model.q_;  // Make backup copy of state
										  // (velocities will be in "v1")

				// k1 = f(t,y);
				// cur_time = t;
				v1 = part->num_model.dotq_;
				// No change needed: part->num_model.q_ = q0;

				part->dyn_simul->can_choose_indep_coords_ = true;
				part->dyn_simul->solve_ddotz(t, ddotz1);
				part->dyn_simul->can_choose_indep_coords_ = false;

				// k2 = f(t+At/2,y+At/2*k1)
				// cur_time = t + t_step2;
				part->dyn_simul->dq_plus_dz(
					v1, t_step2 * ddotz1,
					part->num_model
						.dotq_);  // \dot{q}= \dot{q}_0 + At/2 * \ddot{q}_1
				part->num_model.q_ = q0 + t_step2 * v1;
				part->dyn_simul->correct_dependent_q_dq();

				v2 = part->num_model.dotq_;
				part->dyn_simul->solve_ddotz(t + t_step2, ddotz2);

				// k3 = f(t+At/2,y+At/2*k2)
				// cur_time = t + t_step2;
				// part->num_model.dotq_ = v1 + t_step2*ddotq2;  // \dot{q}=
				// \dot{q}_0 + At/2 * \ddot{q}_2
				part->dyn_simul->dq_plus_dz(
					v1, t_step2 * ddotz2, part->num_model.dotq_);

				part->num_model.q_ = q0 + t_step2 * v2;
				part->dyn_simul->correct_dependent_q_dq();

				v3 = part->num_model.dotq_;
				part->dyn_simul->solve_ddotz(t + t_step2, ddotz3);

				// k4 = f(t+At  ,y+At*k3)
				// cur_time = t + t_step;
				// part->num_model.dotq_ = v1 + t_step*ddotq3;
				part->dyn_simul->dq_plus_dz(
					v1, t_step * ddotz3, part->num_model.dotq_);
				part->num_model.q_ = q0 + t_step * v3;
				part->dyn_simul->correct_dependent_q_dq();

				v4 = part->num_model.dotq_;
				part->dyn_simul->solve_ddotz(t + t_step, ddotz4);

				// Runge-Kutta 4th order formula:
				q_incr = t_step6 * (v1 + 2 * v2 + 2 * v3 + v4);
				dotz_incr =
					t_step6 * (ddotz1 + 2 * ddotz2 + 2 * ddotz3 + ddotz4);
			}

			// generate noise:
			if (dotz_incr.size() != dotz_noise.size())
				dotz_noise.resize(dotz_incr.size());
			random_generator.drawGaussian1DMatrix(
				dotz_noise, 0, model_options.acc_xy_noise_std * t_step);

			// Add (noisy) increment:
			part->num_model.q_ = q0 + q_incr;
			part->dyn_simul->dq_plus_dz(
				v1, dotz_incr + dotz_noise, part->num_model.dotq_);
			part->dyn_simul->correct_dependent_q_dq();

		}  // end numeric integration of one time_step

	}  // end for each time_step

	timelog().leave("PF.1.forward_model");

	// 2) Update weights with sensor measurements:
	// -----------------------------------------------------
	timelog().enter("PF.2.sensor_likelihood");

	const size_t nSensors = sensor_descriptions.size();

	//	std::vector<double> sensors_logw, sensors_loglik;
	//	sensors_logw.reserve( nParts );
	//	sensors_loglik.reserve( nParts );

	for (CParticleList::iterator it = m_particles.begin();
		 it != m_particles.end(); ++it)
	{
		auto part = it->d;

		double cum_log_lik = 0;
		for (size_t k = 0; k < nSensors; k++)
		{
			const double log_lik =
				sensor_descriptions[k]->evaluate_log_likelihood(
					sensor_readings[k], part->num_model);
			cum_log_lik += log_lik;
		}
		it->log_w += cum_log_lik;

		//		sensors_logw.push_back( it->log_w );
		//		sensors_loglik.push_back( cum_log_lik );
	}

	//	double sensor_avrg_lik = mrpt::math::chi2
	// CDF(nSensors,-mrpt::math::averageLogLikelihood(sensors_logw,sensors_loglik)
	//);
	timelog().leave("PF.2.sensor_likelihood");

	//	cout << "Sensor lik: " << sensor_avrg_lik << endl;

	// 3) Normalize weights:
	// ---------------------------------------------------
	timelog().enter("PF.3.renormalize_w");
	this->normalizeWeights();
	timelog().leave("PF.3.renormalize_w");

	// 4) Resampling:
	// -----------------------------------------------------
	timelog().enter("PF.4.resampling");

	const double curESS = this->ESS();
	out_info.resampling_done = false;
	out_info.ESS = curESS;

	if (curESS < PF_options.BETA)
	{
		// printf("[PF] Resampling particles (ESS was %.02f)\n", curESS);

		size_t nNewParts = nParts;	// std::max(40.0, nParts*0.95 );

		this->performResampling(PF_options, nNewParts);	 // Resample

		out_info.resampling_done = true;
	}
	timelog().leave("PF.4.resampling");
}

MultiBodyParticleFilter::TTransitionModelOptions::TTransitionModelOptions()
	: acc_xy_noise_std(1e-3)
{
}

void MultiBodyParticleFilter::getAs3DRepresentation(
	mrpt::opengl::CSetOfObjects::Ptr& outObj,
	const Body::TRenderParams& rp) const
{
	ASSERT_(outObj);

	outObj->clear();
	for (CParticleList::const_iterator it = m_particles.begin();
		 it != m_particles.end(); ++it)
	{
		const auto part = it->d;

		mrpt::opengl::CSetOfObjects::Ptr gl_part =
			mrpt::opengl::CSetOfObjects::Create();
		part->num_model.getAs3DRepresentation(gl_part, rp);
		outObj->insert(gl_part);
	}
}

void MultiBodyParticleFilter::update3DRepresentation(
	const Body::TRenderParams& rp_) const
{
	Body::TRenderParams rp = rp_;

	rp.render_style = Body::reLine;
	for (CParticleList::const_iterator it = m_particles.begin();
		 it != m_particles.end(); ++it)
	{
		const auto part = it->d;

		const uint8_t new_alpha =
			uint8_t(std::max(0.2, std::exp(it->log_w)) * 255);

		rp.line_alpha = new_alpha;

		part->num_model.update3DRepresentation(rp);
	}
}
