/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2024 University of Almeria                           |
  | Copyright (C) 2021 University of Salento                                |
  | See README for list of authors and papers                               |
  | Distributed under 3-clause BSD license                                  |
  |  See: <https://opensource.org/licenses/BSD-3-Clause>                    |
  +-------------------------------------------------------------------------+ */

// Example of particle-filter for tracking a multibody model
// ------------------------------------------------------------
#include <mbse/mbse.h>
#include <mbse/MultiBodyParticleFilter.h>
#include <mbse/constraints/ConstraintMobileSlider.h>
#include <mbse/constraints/ConstraintFixedSlider.h>

#include <mrpt/opengl.h>
#include <mrpt/gui.h>
#include <mrpt/poses/CPosePDFParticles.h>
#include <mrpt/io/vector_loadsave.h>

#include <thread>  // for sleep()

using namespace std;
using namespace mbse;
using namespace mrpt;
using namespace mrpt::poses;
using namespace mrpt::math;

// ===================================================================
//						EXPERIMENT CONFIGURATION
// ===================================================================
const size_t NUM_PARTS =
	getenv("PF_NUM_SAMPLES") != nullptr ? atoi(getenv("PF_NUM_SAMPLES")) : 500;

// PF's sensor model: 1 sigma of the assumed noise (in meters)
const double PF_MOTION_MODEL_NOISE_XY =
	getenv("PF_MOTION_MODEL_NOISE_XY") != nullptr
		? atof(getenv("PF_MOTION_MODEL_NOISE_XY"))
		: 1;

// Gyro's assumed noise (1 sigma) (deg/s)
const double PF_SENSOR_NOISE_STD = DEG2RAD(
	getenv("PF_SENSOR_NOISE_STD") != nullptr
		? atof(getenv("PF_SENSOR_NOISE_STD"))
		: 10.0);  // PF's sensor model: 1 sigma of the assumed noise (in meters)

// Real noise added to simulated sensor
const double REAL_SENSOR_STD = DEG2RAD(0.3);

// Time step lenght (in seconds):
const double SIMUL_STEP =
	getenv("SIMUL_STEP") != nullptr ? atof(getenv("SIMUL_STEP")) : 5e-3;

const int DRAW_DECIMATION = 1;
const int DRAW_DELAY_MS =
	getenv("DRAW_DELAY_MS") != nullptr ? atoi(getenv("DRAW_DELAY_MS")) : 1;

// Optionally: modify mass of body #1 by a factor different than 1.0 to force
// the PF to track an incorrectly modeled system:
const double IMPERFECT_PF_MODEL_ERROR =
	getenv("IMPERFECT_PF_MODEL_ERROR") != nullptr
		? atof(getenv("IMPERFECT_PF_MODEL_ERROR"))
		: 1;

const bool INITIALIZE_UNIFORMLY = true;

const double FINAL_TIME = 10.0;

#define SHOW_GUI
#define COMPUTE_RMSE 1
#define SAVE_STATS 0

// ===================================================================
// ===================================================================

void buildFourBarsMBS(ModelDefinition& model)
{
	Body::TRenderParams rp;
	rp.render_style = Body::reLine;
	rp.line_width = 3.0f;
	// rp.cyl_diameter = 0.01;

	model.setPointCount(4);
	model.setPointCoords(0, TPoint2D(0, 0), true /*is fixed*/);
	model.setPointCoords(1, TPoint2D(0.125, 0));
	model.setPointCoords(2, TPoint2D(0.125, 0.270));
	model.setPointCoords(3, TPoint2D(0.590, 0), true /*is fixed*/);

	{
		Body& b = model.addBody();
		b.points[0] = 0;
		b.points[1] = 1;

		b.length() = 0.125;
		b.mass() = 39e-3;
		b.I0() = (1. / 3.) * b.mass() * square(b.length());
		b.cog() = TPoint2D(b.length() * 0.5, 0);
		b.render_params = rp;
	}
	{
		Body& b = model.addBody();
		b.points[0] = 1;
		b.points[1] = 2;

		b.length() = 0.270;
		b.mass() = 161.13e-3;
		b.I0() = (1. / 3.) * b.mass() * square(b.length());
		b.cog() = TPoint2D(b.length() * 0.5, 0);
		b.render_params = rp;
	}
	{
		Body& b = model.addBody();
		b.points[0] = 2;
		b.points[1] = 3;

		b.length() = 0.540;	 // std::sqrt(2.0*2.0+3.0*3.0);
		b.mass() = 303.70e-3;
		b.I0() = (1. / 3.) * b.mass() * square(b.length());
		b.cog() = TPoint2D(b.length() * 0.5, 0);
		b.render_params = rp;
	}
}

void buildFollowerMBS(ModelDefinition& model)
{
	model.setPointCount(5);
	model.setPointCoords(0, TPoint2D(0, 0), true /*is fixed*/);
	model.setPointCoords(1, TPoint2D(1, 1));
	model.setPointCoords(2, TPoint2D(1, -4), true /*is fixed*/);
	model.setPointCoords(3, TPoint2D(1, 4));
	model.setPointCoords(4, TPoint2D(5, 0));

	{
		Body& b = model.addBody();
		b.points[0] = 0;
		b.points[1] = 1;

		b.length() = std::sqrt(2);
		b.mass() = 1;
		b.I0() = (1. / 3.) * b.mass() * square(b.length());
		b.cog() = TPoint2D(b.length() * 0.5, 0);
	}
	{
		Body& b = model.addBody();
		b.points[0] = 2;
		b.points[1] = 3;

		b.length() = 8;
		b.mass() = 3;
		b.I0() = (1. / 3.) * b.mass() * square(b.length());
		b.cog() = TPoint2D(b.length() * 0.5, 0);
	}
	{
		Body& b = model.addBody();
		b.points[0] = 3;
		b.points[1] = 4;

		b.length() = std::sqrt(32);
		b.mass() = 2;
		b.I0() = (1. / 3.) * b.mass() * square(b.length());
		b.cog() = TPoint2D(b.length() * 0.5, 0);
	}

	model.addConstraint<ConstraintMobileSlider>(
		1 /*pt index*/, 2, 3 /* two pts for defining the constraint */
	);

	model.addConstraint<ConstraintFixedSlider>(
		4 /*pt index*/, TPoint2D(-5, 0),
		TPoint2D(10, 0) /* The line on which to fix the point */
	);
}

void pf_initialize_uniform_distribution(
	MultiBodyParticleFilter& pf, ModelDefinition& model);

auto& rnd = mrpt::random::getRandomGenerator();

int main(int argc, char** argv)
{
	try
	{
		mbse::timelog().enable(false);

		ModelDefinition model;

		// Define the mechanism:
		// ------------------------------------------------
		buildFourBarsMBS(model);
		// buildFollowerMBS(model);

		// "Compile" the problem, and use as "Ground Truth" (GT) for the
		// estimation problem
		// -------------------------------------------------------------------------------------
		std::shared_ptr<AssembledRigidModel> aMBS_GT = model.assembleRigidMBS();

		// Prepare ground truth dynamic simulation:
		// -----------------------------------------------
		CDynamicSimulator_Lagrange_LU_dense dynSimul_GT(aMBS_GT);

		// Set params:
		// -----------------------------
		dynSimul_GT.params.time_step = SIMUL_STEP;
		dynSimul_GT.params.ode_solver = ODE_RK4;

		// Prepare solver; must be called before "run()".
		// ----------------------------------------------
		dynSimul_GT.prepare();

		// Prepare Particle Filter:
		// -----------------------------------------------
		ModelDefinition PF_model;
		buildFourBarsMBS(PF_model);
		// buildFollowerMBS(PF_model);

		// Corrupt PF model:
		{
			const double error = IMPERFECT_PF_MODEL_ERROR;

			std::vector<Body>& bodies = PF_model.bodies();
			bodies[1].length() *= error;
			bodies[1].mass() *= error;
			bodies[1].I0() =
				(1. / 3.) * bodies[1].mass() * square(bodies[1].length());
			bodies[1].cog() = TPoint2D(bodies[1].length() * 0.5, 0);
		}

		MultiBodyParticleFilter pf(NUM_PARTS, PF_model);

		// PF parameters:
		// ------------------------
		pf.model_options.acc_xy_noise_std = PF_MOTION_MODEL_NOISE_XY;
		// pf.PF_options.verbose = true;

		// PF must be initialized with a uniform distribution (no "prior" at
		// all)
		//  For now, it's done BY HAND for this specific mechanism...
		// -------------------------------------------------------------------------

		if (INITIALIZE_UNIFORMLY)
		{
			mbse::timelog().enter("PF_unif_initialize");
			pf_initialize_uniform_distribution(pf, model);
			mbse::timelog().leave("PF_unif_initialize");
		}

#ifdef SHOW_GUI
		// Prepare 3D scene:
		// -----------------------------------------------
		Body::TRenderParams rp;
		Body::TRenderParams rp_particles;
		rp_particles.render_style = Body::reLine;

		mrpt::opengl::CSetOfObjects::Ptr gl_MBS =
			mrpt::opengl::CSetOfObjects::Create();
		aMBS_GT->getAs3DRepresentation(gl_MBS, rp);

		mrpt::opengl::CSetOfObjects::Ptr gl_MBPF =
			mrpt::opengl::CSetOfObjects::Create();
		pf.getAs3DRepresentation(gl_MBPF, rp_particles);

		mrpt::gui::CDisplayWindow3D win3D("MBS dynamic simulation", 1280, 720);

		win3D.setCameraAzimuthDeg(-90);
		win3D.setCameraElevationDeg(90);
		win3D.setCameraZoom(1);
		mrpt::opengl::COpenGLViewport::Ptr gl_estimate;
		{
			auto& scene = win3D.get3DSceneAndLock();

			auto gl_main = scene->getViewport("main");
			gl_estimate = scene->createViewport("estimation");

			// Main view:
			gl_main->setViewportPosition(0, 0.5, 1.0, 0.5);
			gl_main->setBorderSize(1);
			gl_main->insert(
				mrpt::opengl::stock_objects::CornerXYZSimple(0.2, 1));
			gl_main->insert(gl_MBS);

			mrpt::opengl::TFontParams fp;
			fp.color = mrpt::img::TColorf(0, 0, 1);
			fp.vfont_scale = 20;  // pixels

			win3D.addTextMessage(10, 0.95, "Ground truth", 1 /* txt ID */, fp);

			// Estimate:
			gl_estimate->setViewportPosition(0, 0, 1.0, 0.5);
			gl_estimate->setBorderSize(1);
			gl_estimate->insert(
				mrpt::opengl::stock_objects::CornerXYZSimple(0.2, 1));
			gl_estimate->insert(gl_MBPF);

			win3D.addTextMessage(10, 0.45, "Estimation", 2 /* txt ID */, fp);

			win3D.unlockAccess3DScene();
		}
		win3D.repaint();

		// Save initial state:
#if 0
		{
			mrpt::opengl::COpenGLScene scene;
			scene.insert(gl_MBPF);
			scene.saveToFile("initial.3Dscene");
		}
#endif

		// win3D.grabImagesStart();
#endif	// SHOW_GUI

#if SAVE_STATS || COMPUTE_RMSE
		vector<double> STATS_t, GT_ang0, EST_ang0_mean, EST_ang0_std;
		mrpt::poses::CPosePDFParticles orientation_averager;
		orientation_averager.resetDeterministic(
			mrpt::math::TPose2D(), pf.m_particles.size());
#endif

		// Prepare sensor descriptions:
		// ---------------------------------------
		std::vector<CVirtualSensor::Ptr> sensor_descriptions;
		std::vector<double> sensor_readings;

		// Gyroscope on body #2:
		//		sensor_descriptions.push_back( CVirtualSensor::Ptr(new
		// CVirtualSensor_Gyro(2 /* body index */)) );
		//		sensor_descriptions.back()->sensor_noise_std =
		// PF_SENSOR_NOISE_STD;

		// Gyroscope on body #1:
		sensor_descriptions.push_back(
			CVirtualSensor::Ptr(new CVirtualSensor_Gyro(1 /* body index */)));
		sensor_descriptions.back()->sensor_noise_std = PF_SENSOR_NOISE_STD;

		// -----------------------------------------------
		//         Run PF estimator experiment
		// -----------------------------------------------
		double t_old_simul = 0;

		while (t_old_simul < FINAL_TIME)
		{
			// Run the "real" (ground truth) model:
			// ==============================================
			const double t_pf_start = t_old_simul;
			{
				const double t_new_simul = t_old_simul + SIMUL_STEP;
				// Simulate: (Return the actual end-time of the simulation,
				// >=t_new)
				const double t_final =
					dynSimul_GT.run(t_old_simul, t_new_simul);
				t_old_simul = t_final;
			}
			const double t_pf_end = t_old_simul;

			// Simulate "real" sensors:
			const size_t nSensors = sensor_descriptions.size();
			sensor_readings.resize(nSensors);
			for (size_t j = 0; j < nSensors; j++)
				sensor_readings[j] =
					sensor_descriptions[j]->simulate_reading(*aMBS_GT) +
					rnd.drawGaussian1D(0, REAL_SENSOR_STD);

			// Run the PF models:
			// ==============================================
			MultiBodyParticleFilter::TOutputInfo pf_out_info;
			pf.run_PF_step(
				t_pf_start, t_pf_end, SIMUL_STEP, sensor_descriptions,
				sensor_readings, pf_out_info);

#if SAVE_STATS || COMPUTE_RMSE
			// Save stats of the estimation quality:
			// =====================================================
			STATS_t.push_back(t_old_simul);
			// GT value:
			GT_ang0.push_back(atan2(
				aMBS_GT->q_[aMBS_GT->points2DOFs_[1].dof_y],
				aMBS_GT->q_[aMBS_GT->points2DOFs_[1].dof_x]));

			// Estimated values:
			for (size_t i = 0; i < pf.m_particles.size(); i++)
			{
				auto part = pf.m_particles[i].d;

				const double val_phi = atan2(
					part->num_model.q_[part->num_model.points2DOFs_[1].dof_y],
					part->num_model.q_[part->num_model.points2DOFs_[1].dof_x]);

				orientation_averager.m_particles[i].d.phi = val_phi;

				orientation_averager.m_particles[i].log_w =
					pf.m_particles[i].log_w;
			}

			// Estimate a full (x,y,phi) pose just for easy reusing of MRPT
			// implementation:
			const auto [parts_cov, parts_mean] =
				orientation_averager.getCovarianceAndMean();

			// and only keep the orientation:
			EST_ang0_mean.push_back(parts_mean.phi());
			EST_ang0_std.push_back(std::sqrt(parts_cov(2, 2)));
#endif

			// After first resampling, set actual sensor noise:
			if (pf_out_info.resampling_done)
			{
				static int FACTOR = 10;
				if (FACTOR > 1) --FACTOR;

				for (size_t k = 0; k < sensor_descriptions.size(); k++)
					sensor_descriptions[k]->sensor_noise_std =
						PF_SENSOR_NOISE_STD * FACTOR;
			}

#ifdef SHOW_GUI
			// Upon resampling we have to rebuild the 3D rendering:
			if (pf_out_info.resampling_done)
			{
				win3D.get3DSceneAndLock();

				gl_MBPF->clear();
				pf.getAs3DRepresentation(gl_MBPF, rp_particles);

				win3D.unlockAccess3DScene();
			}

			// Handle key-strokes:
			if (win3D.keyHit())
			{
				const int c = win3D.getPushedKey();
				switch (c)
				{
						// Re-initialize:
					case 'R':
					case 'r':
						mbse::timelog().enter("PF_unif_initialize");
						pf_initialize_uniform_distribution(pf, model);
						mbse::timelog().leave("PF_unif_initialize");
						break;
				};
			}

			// Update 3D view:
			// ==============================================
			{
				static int draw_decim_cnt = 0;
				if (++draw_decim_cnt >= DRAW_DECIMATION)
				{
					draw_decim_cnt = 0;

					// DRAW:
					win3D.get3DSceneAndLock();

					aMBS_GT->update3DRepresentation(rp);
					pf.update3DRepresentation(rp_particles);

					// Replicate camera view in both viewports:
					mrpt::opengl::CCamera& cam = gl_estimate->getCamera();
					cam.setAzimuthDegrees(win3D.getCameraAzimuthDeg());
					cam.setElevationDegrees(win3D.getCameraElevationDeg());
					cam.setZoomDistance(win3D.getCameraZoom());
					float px, py, pz;
					win3D.getCameraPointingToPoint(px, py, pz);
					cam.setPointingAt(px, py, pz);

					win3D.unlockAccess3DScene();
					win3D.repaint();
				}

				// Update 3D scene:

				mrpt::opengl::TFontParams fp;
				fp.color = mrpt::img::TColorf(1, 1, 1);
				fp.vfont_scale = 15;  // pixels

				win3D.addTextMessage(
					10, 10,
					mrpt::format(
						"Time: %.03fs | ESS=%.02f | %u samples", t_old_simul,
						pf_out_info.ESS,
						static_cast<unsigned int>(pf.particlesCount())),
					0 /* txt ID */, fp);

				std::this_thread::sleep_for(
					std::chrono::milliseconds(DRAW_DELAY_MS));
			}
#endif	// SHOW_GUI
		}

#if SAVE_STATS
		mrpt::io::vectorToTextFile(STATS_t, "t.txt");
		mrpt::io::vectorToTextFile(GT_ang0, "GT_ang0.txt");
		mrpt::io::vectorToTextFile(EST_ang0_mean, "EST_ang0_mean.txt");
		mrpt::io::vectorToTextFile(EST_ang0_std, "EST_ang0_std.txt");
#endif

#if COMPUTE_RMSE
		{
			double sqerr = 0;
			size_t n = 0;
			for (size_t i = STATS_t.size() * 0.75; i < STATS_t.size(); i++)
			{
				sqerr += square(
					mrpt::math::wrapTo2Pi(GT_ang0[i]) -
					mrpt::math::wrapTo2Pi(EST_ang0_mean[i]));
				n++;
			}
			const double rmse = std::sqrt(sqerr / n);
			cout << rmse << endl;
		}
#endif

		return 0;  // program ended OK.
	}
	catch (exception& e)
	{
		cerr << e.what() << endl;
		return 1;
	}
}

void pf_initialize_uniform_distribution(
	MultiBodyParticleFilter& pf, ModelDefinition& model)
{
	for (size_t i = 0; i < pf.m_particles.size(); i++)
	{
		auto part = pf.m_particles[i].d;
		const auto& p2dofs = part->num_model.getPoints2DOFs();

		// We draw a random value for the free DOF(s) of the mechanism:
		// and uniform values for the rest:
		double final_err;
		do
		{
			const double R = model.bodies()[0].length();
			const double ang = rnd.drawUniform(-M_PI, M_PI);
			const size_t nTotalDOFs = part->num_model.q_.size();
			for (size_t k = 0; k < nTotalDOFs; k++)
				part->num_model.q_[k] = rnd.drawUniform(-20, 20);

			const size_t PT_IDX = 1;  // This point is the one we force to be at
									  // a predefined position:
			part->num_model.q_[p2dofs[PT_IDX].dof_x] = R * cos(ang);
			part->num_model.q_[p2dofs[PT_IDX].dof_y] = R * sin(ang);

			final_err = part->num_model.refinePosition(1e-13, 30);
		} while (final_err > 1e-6);
	}
}
