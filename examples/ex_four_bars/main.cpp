
// Example of MBS model: Double Quad
// -------------------------------------
#include <sparsembs/sparsembs.h>
#include <sparsembs/model-examples.h>

#include <mrpt/opengl.h>
#include <mrpt/gui.h>
#include <mrpt/math/ops_vectors.h>
#include <thread>  // for sleep()

using namespace std;
using namespace sparsembs;
using namespace mrpt;
using namespace mrpt::poses;
using namespace mrpt::math;

#define SIMUL_REALTIME 1  // 1: real time, 0: fixed simulation step

double REALTIME_FACTOR = 0.15;

const double GUI_DESIRED_FPS = 75;	// Hz

const bool SHOW_ENERGY_BALANCE = true;

void my_callback(const TSimulationStateRef simul_state) {}

<<<<<<< HEAD
=======
void buildParameterizedMBS(
	const size_t nx, const size_t ny, CModelDefinition& model,
	const double NOISE_LEN = 0)
{
	ASSERT_(nx >= 1 && ny >= 1);

#if MRPT_VERSION >= 0x199
	auto randomGenerator = mrpt::random::getRandomGenerator();
#else
	auto& randomGenerator = mrpt::random::randomGenerator;
#endif

	// Definition of constants
	const size_t npoints = (nx + 1) * (ny + 1);  // number of points

	const double Lx = 2.0;  // Rod lengths
	const double Ly = 2.5;  // Rod lengths

	model.setPointCount(npoints);

	// Definition of fixed points
	for (size_t i = 0; i <= nx; i++)
	{
		const double x =
			i * Lx;  // + (irregular_mesh ?
					 // randomGenerator.drawUniform(-NOISE_LEN,NOISE_LEN) : 0 );
		model.setPointCoords(i, TPoint2D(x, 0), true /*is fixed*/);
	}

	// points definition
	for (size_t row = 1; row <= ny;
		 row++)  // from row 1 (0 belongs to fixed points) to ny
	{
		for (size_t j = 0; j <= nx; j++)  // from column 0 to nx
		{
			const double x =
				j * Lx + randomGenerator.drawUniform(-NOISE_LEN, NOISE_LEN);
			const double y =
				row * Ly + randomGenerator.drawUniform(-NOISE_LEN, NOISE_LEN);
			model.setPointCoords(row * (nx + 1) + j, TPoint2D(x, y));
		}
	}

	// horizontal bars
	for (size_t row = 1; row <= ny;
		 row++)  // from row 1 (0 belongs to fixed points) to ny
	{
		for (size_t j = 0; j < nx; j++)  // from column 0 to nx-1
		{
			CBody& b = model.addBody();
			b.points[0] = row * (nx + 1) + j;  // left side point
			b.points[1] = row * (nx + 1) + j + 1;  // right side point
			b.length() = (model.getPointInfo(b.points[0]).coords -
						  model.getPointInfo(b.points[1]).coords)
							 .norm();
			b.mass() = 1 * b.length();
			b.I0() = b.mass() * mrpt::square(b.length()) / 3.0;
			b.cog() = TPoint2D(b.length() * 0.5, 0);
		}
	}

	// vertical bars
	for (size_t row = 1; row <= ny;
		 row++)  // from row 1 (0 belongs to fixed points) to ny
	{
		for (size_t j = 0; j <= nx; j++)  // from column 0 to nx
		{
			CBody& b = model.addBody();
			b.points[0] = row * (nx + 1) + j;  // upper point
			b.points[1] = row * (nx + 1) + j - nx - 1;  // lower point
			b.length() = (model.getPointInfo(b.points[0]).coords -
						  model.getPointInfo(b.points[1]).coords)
							 .norm();
			b.mass() = 1 * b.length();
			b.I0() = b.mass() * mrpt::square(b.length()) / 3.0;
			b.cog() = TPoint2D(b.length() * 0.5, 0);

			b.render_params.z_layer = 0.1;
		}
	}

	cout << "# bodies: " << model.getBodies().size() << endl;
}

void buildLongStringMBS(const size_t N, CModelDefinition& model)
{
	ASSERT_(N >= 1);

	// Definition of constants
	const double L = 0.5;  // Rod lengths

	model.setPointCount(N + 1);

	// Definition of fixed points
	for (size_t i = 0; i <= N; i++)
		model.setPointCoords(i, TPoint2D(i * L, 0), i == 0 /*is fixed*/);

	// bars
	for (size_t j = 0; j < N; j++)
	{
		CBody& b = model.addBody();
		b.points[0] = j;
		b.points[1] = j + 1;
		b.length() = L;
		b.mass() = L * 0.1;
		b.I0() = b.mass() * mrpt::square(b.length()) / 3.0;
		b.cog() = TPoint2D(b.length() * 0.5, 0);
	}
}

void buildFourBarsMBS(CModelDefinition& model)
{
	model.setPointCount(4);
	model.setPointCoords(0, TPoint2D(0, 0), true /*is fixed*/);
	model.setPointCoords(1, TPoint2D(1, 0));
	model.setPointCoords(2, TPoint2D(1, 2));
	model.setPointCoords(3, TPoint2D(4, 0), true /*is fixed*/);

	{
		CBody& b = model.addBody();
		b.points[0] = 0;
		b.points[1] = 1;

		b.length() = 1;
		b.mass() = 1;
		b.I0() = (1. / 3.) * b.mass() * square(b.length());
		b.cog() = TPoint2D(b.length() * 0.5, 0);

		b.render_params.z_layer = 0;
	}
	{
		CBody& b = model.addBody();
		b.points[0] = 1;
		b.points[1] = 2;

		b.length() = 2;
		b.mass() = 2;
		b.I0() = (1. / 3.) * b.mass() * square(b.length());
		b.cog() = TPoint2D(b.length() * 0.5, 0);

		b.render_params.z_layer = -0.05;
	}
	{
		CBody& b = model.addBody();
		b.points[0] = 2;
		b.points[1] = 3;

		b.length() = std::sqrt(2.0 * 2.0 + 3.0 * 3.0);
		b.mass() = 4;
		b.I0() = (1. / 3.) * b.mass() * square(b.length());
		b.cog() = TPoint2D(b.length() * 0.5, 0);

		b.render_params.z_layer = 0;
	}
}

void buildFourBarsMBS_JavierCuadrado(CModelDefinition& model)
{
	model.setPointCount(4);
	model.setPointCoords(0, TPoint2D(0, 0), true /*is fixed*/);
	model.setPointCoords(1, TPoint2D(1, 1.732050807568877));
	model.setPointCoords(2, TPoint2D(8.412459326544008, 4.741277740242908));
	model.setPointCoords(3, TPoint2D(10, 0), true /*is fixed*/);

	{
		CBody& b = model.addBody();
		b.points[0] = 0;
		b.points[1] = 1;

		b.length() = 2;
		b.mass() = 2;
		b.I0() = (1. / 3.) * b.mass() * square(b.length());
		b.cog() = TPoint2D(b.length() * 0.5, 0);
	}
	{
		CBody& b = model.addBody();
		b.points[0] = 1;
		b.points[1] = 2;

		b.length() = 8;
		b.mass() = 8;
		b.I0() = (1. / 3.) * b.mass() * square(b.length());
		b.cog() = TPoint2D(b.length() * 0.5, 0);
	}
	{
		CBody& b = model.addBody();
		b.points[0] = 2;
		b.points[1] = 3;

		b.length() = 5;
		b.mass() = 5;
		b.I0() = (1. / 3.) * b.mass() * square(b.length());
		b.cog() = TPoint2D(b.length() * 0.5, 0);
	}
}

void buildSliderCrankMBS(CModelDefinition& model)
{
	model.setPointCount(3);
	model.setPointCoords(0, TPoint2D(0, 0), true /*is fixed*/);
	model.setPointCoords(1, TPoint2D(1, 1));
	model.setPointCoords(2, TPoint2D(5, 0));

	{
		CBody& b = model.addBody();
		b.points[0] = 0;
		b.points[1] = 1;

		b.length() = std::sqrt(2);
		b.mass() = 1;
		b.I0() = (1. / 3.) * b.mass() * square(b.length());
		b.cog() = TPoint2D(b.length() * 0.5, 0);
	}
	{
		CBody& b = model.addBody();
		b.points[0] = 1;
		b.points[1] = 2;

		b.length() = std::sqrt(17);
		b.mass() = 2;
		b.I0() = (1. / 3.) * b.mass() * square(b.length());
		b.cog() = TPoint2D(b.length() * 0.5, 0);
	}

	model.addConstraint(CConstraintFixedSlider(
		2 /*pt index*/, TPoint2D(-3, -2),
		TPoint2D(8, 2) /* The line on which to fix the point */
		));
}

void buildFollowerMBS(CModelDefinition& model)
{
	model.setPointCount(5);
	model.setPointCoords(0, TPoint2D(0, 0), true /*is fixed*/);
	model.setPointCoords(1, TPoint2D(1, 1));
	model.setPointCoords(2, TPoint2D(1, -4), true /*is fixed*/);
	model.setPointCoords(3, TPoint2D(1, 4));
	model.setPointCoords(4, TPoint2D(5, 0));

	{
		CBody& b = model.addBody();
		b.points[0] = 0;
		b.points[1] = 1;

		b.length() = std::sqrt(2);
		b.mass() = 1;
		b.I0() = (1. / 3.) * b.mass() * square(b.length());
		b.cog() = TPoint2D(b.length() * 0.5, 0);
	}
	{
		CBody& b = model.addBody();
		b.points[0] = 2;
		b.points[1] = 3;

		b.length() = 8;
		b.mass() = 3;
		b.I0() = (1. / 3.) * b.mass() * square(b.length());
		b.cog() = TPoint2D(b.length() * 0.5, 0);
	}
	{
		CBody& b = model.addBody();
		b.points[0] = 3;
		b.points[1] = 4;

		b.length() = std::sqrt(32);
		b.mass() = 2;
		b.I0() = (1. / 3.) * b.mass() * square(b.length());
		b.cog() = TPoint2D(b.length() * 0.5, 0);
	}

	model.addConstraint(CConstraintMobileSlider(
		1 /*pt index*/, 2, 3 /* two pts for defining the constraint */
		));

	model.addConstraint(CConstraintFixedSlider(
		4 /*pt index*/, TPoint2D(-5, 0),
		TPoint2D(10, 0) /* The line on which to fix the point */
		));
}

void buildTwoSliderBlocks(CModelDefinition& model)
{
	model.setPointCount(2);
	model.setPointCoords(0, TPoint2D(0, 15 * sin(DEG2RAD(35))));
	model.setPointCoords(1, TPoint2D(15 * cos(DEG2RAD(35)), 0));

	{
		CBody& b = model.addBody();
		b.points[0] = 0;
		b.points[1] = 1;

		b.length() = 15;
		b.mass() = 1;
		b.I0() = (1. / 3.) * b.mass() * square(b.length());
		b.cog() = TPoint2D(b.length() * 0.5, 0);
	}

	model.addConstraint(CConstraintFixedSlider(
		0 /*pt index*/, TPoint2D(0, 0),
		TPoint2D(0, 1) /* The line on which to fix the point */
		));
	model.addConstraint(CConstraintFixedSlider(
		1 /*pt index*/, TPoint2D(0, 0),
		TPoint2D(1, 0) /* The line on which to fix the point */
		));
}

>>>>>>> origin/master
int main(int argc, char** argv)
{
	try
	{
		CBody::TRenderParams dynamic_rp;

		CModelDefinition model;

		// Define a double pendulum mechanism:
		// ------------------------------------------------
		// buildFourBarsMBS(model);
		// buildFourBarsMBS_JavierCuadrado(model);
		// buildSliderCrankMBS(model);
		// buildFollowerMBS(model);
		// buildLongStringMBS( 15, model);
		sparsembs::buildTwoSliderBlocks(model);

		// const size_t Nx = 5, Ny = 4;
		// buildParameterizedMBS(Nx, Ny, model, 0.6 /*random imperfections in
		// the mesh*/);

		// "Compile" the problem:
		// ------------------------------------------------
		std::shared_ptr<CAssembledRigidModel> aMBS = model.assembleRigidMBS();

		// Set initial velocities (Only for the buildParameterizedMBS model)
		// for (size_t i=0;i<=Nx;i++)
		//	aMBS->m_dotq[2*i+0]=1;

#if 0
		// prueba del ejercicio de clase (model: buildTwoSliderBlocks)
		{
			aMBS->m_dotq[1]=-10;

			std::vector<size_t> z_indices;
			Eigen::VectorXd out_ddq, ddot_z;
			z_indices.push_back(1);
			CAssembledRigidModel::TComputeDependentParams cdp;
			CAssembledRigidModel::TComputeDependentResults res;
			out_ddq.resize(4);
			res.ddotq = &out_ddq;

			ddot_z.resize(1);
			ddot_z[0] = -5;

			aMBS->computeDependentPosVelAcc(z_indices,true,true,cdp,res,&ddot_z);
		}
#endif

		// Change gravity vector value:
		// aMBS->setGravityVector(0,-9.80665,0);
		aMBS->setGravityVector(0, -9.81, 0);

		// Prepare 3D scene:
		// -----------------------------------------------
		mrpt::opengl::CSetOfObjects::Ptr gl_MBS =
			mrpt::opengl::CSetOfObjects::Create();
		aMBS->getAs3DRepresentation(gl_MBS, dynamic_rp);

		mrpt::gui::CDisplayWindow3D win3D("MBS dynamic simulation", 1280, 768);
		win3D.setCameraAzimuthDeg(90);
		win3D.setCameraElevationDeg(0);
		win3D.setCameraZoom(45);

		// win3D.setCameraPointingToPoint(Nx*2*0.5,0, Ny*2.5*0.5);

		{
			auto& scene = win3D.get3DSceneAndLock();
			// scene->getViewport()->setCustomBackgroundColor(
			// mrpt::utils::TColorf(0.1f,0.1f,0.1f));

			gl_MBS->setPose(mrpt::poses::CPose3D(
				0, 0, 0, DEG2RAD(180), DEG2RAD(0), DEG2RAD(90)));

			// scene->insert( mrpt::opengl::stock_objects::CornerXYZ() );
			scene->insert(gl_MBS);

			win3D.unlockAccess3DScene();
		}
		win3D.repaint();

		const double GUI_DESIRED_PERIOD = 1.0 / GUI_DESIRED_FPS;
		CTicTac tictac_gui_refresh;

		// Test initial position problem:
		// ---------------------------------
#if 0
		{
			cout << "Press any key to solve initial pos problem\n";
			win3D.waitForKey();
			// Solve and display again:
			const double final_Phi = aMBS->solveInitialPosition();
			cout << "Final |Phi| = " << final_Phi << endl;
			{
				win3D.get3DSceneAndLock();
				aMBS->update3DRepresentation();
				win3D.unlockAccess3DScene();
				win3D.repaint();
			}
			win3D.waitForKey();
		}
#endif

		// Executes a dynamic simulation:
		// -----------------------------------------------
		// CDynamicSimulator_Lagrange_LU_dense
		// CDynamicSimulator_Lagrange_UMFPACK
		// CDynamicSimulator_Lagrange_KLU
		// CDynamicSimulator_Lagrange_CHOLMOD
		//CDynamicSimulator_Indep_dense
		// CDynamicSimulator_AugmentedLagrangian_KLU
		// CDynamicSimulator_AugmentedLagrangian_Dense
		// CDynamicSimulator_ALi3_Dense

		CDynamicSimulator_R_matrix_dense dynSimul(aMBS);

		// dynSimul.params_penalty.alpha = 1e7;
		// dynSimul.params_penalty.xi = 1;
		// dynSimul.params_penalty.w = 10;

		// Mark points for logging:
		// -----------------------------------------------
		// dynSimul.addPointSensor(3);
		// dynSimul.addPointSensor(2);

		// Set params:
		// -----------------------------
		dynSimul.params.time_step = 1e-3;
		dynSimul.params.ode_solver = ODE_RK4;  // ODE_Trapezoidal; //ODE_RK4; //
		dynSimul.params.user_callback = simul_callback_t(my_callback);

		// Energy stats:
		CAssembledRigidModel::TEnergyValues energy;
		vector<double> E_tot, E_kin, E_pot;
		const size_t ENERGY_MAX_LOG = 100000;  // Not to slow down graphs...
		E_tot.reserve(ENERGY_MAX_LOG);
		E_kin.reserve(ENERGY_MAX_LOG);
		E_pot.reserve(ENERGY_MAX_LOG);

		const size_t ENERGY_DECIMATE = 100;
		size_t ENERGY_DECIMATE_CNT = 0;

		// Prepare solver; must be called before "run()".
		// ----------------------------------------------
		dynSimul.prepare();

		// Run dynamic simulation:
		CTicTac tictac;
		double t_old = tictac.Tac();
		double t_old_simul = t_old;

		while (win3D.isOpen())
		{
			double t_new = tictac.Tac();

#if SIMUL_REALTIME
			if (t_new >
				t_old + dynSimul.params.time_step)	// Just in case the computer
													// is *really fast*...
			{
				double t_new_simul =
					t_old_simul + REALTIME_FACTOR * (t_new - t_old);
				// Simulate: (Return the actual end-time of the simulation,
				// >=t_new)
				t_old_simul = dynSimul.run(t_old_simul, t_new_simul);
				t_old = t_new;
			}
#else
			{
				double t_new_simul = t_old_simul + dynSimul.params.time_step;
				t_old_simul = dynSimul.run(t_old_simul, t_new_simul);
				t_old = t_new;
			}
#endif

			// Eval energy:
			if (++ENERGY_DECIMATE_CNT > ENERGY_DECIMATE &&
				E_tot.size() < ENERGY_MAX_LOG)
			{
				ENERGY_DECIMATE_CNT = 0;

				aMBS->evaluateEnergy(energy);
				E_tot.push_back(energy.E_total);
				E_kin.push_back(energy.E_kin);
				E_pot.push_back(energy.E_pot);
			}

			// Update 3D view:
			if (tictac_gui_refresh.Tac() >= GUI_DESIRED_PERIOD)
			{
				tictac_gui_refresh.Tic();

				sparsembs::timelog.enter("update_3D_view");
				win3D.get3DSceneAndLock();
				aMBS->update3DRepresentation(dynamic_rp);
				win3D.unlockAccess3DScene();
				win3D.repaint();
				sparsembs::timelog.leave("update_3D_view");

				// Update 3D scene:
				win3D.addTextMessage(
					10, 10,
					mrpt::format(
						"Time: %.03fs (x%.02f) |Phi|=%e", t_old,
						REALTIME_FACTOR, aMBS->m_Phi.norm()),
					TColorf(1, 1, 1), "mono", 15, mrpt::opengl::NICE,
					0 /* txt ID */, 1.5 /*spacing*/, 0.1 /*kerning*/,
					true /*draw shadow */);

				const double simul_t =
					sparsembs::timelog.getMeanTime("mbs.run_complete_timestep");
				const double simul_Hz = simul_t > 0 ? 1.0 / simul_t : 0;
				win3D.addTextMessage(
					10, 30,
					mrpt::format(
						"Simul: %.02fHz [%.02fms]", simul_Hz, simul_t * 1e3),
					TColorf(1, 1, 1), "mono", 15, mrpt::opengl::NICE,
					1 /* txt ID */, 1.5 /*spacing*/, 0.1 /*kerning*/,
					true /*draw shadow */);

				if (win3D.keyHit())
				{
					const int k = win3D.getPushedKey();

					switch (k)
					{
#if SIMUL_REALTIME
						case '+':
							REALTIME_FACTOR *= 1.1;
							break;
						case '-':
							REALTIME_FACTOR /= 1.1;
							break;
#endif
						default:
							break;
					};
				}

			}  // end update GUI

			std::this_thread::sleep_for(std::chrono::milliseconds(10));
		}

		// Save sensed data:
		// ----------------------------
		cout << "Saving log data to LOG.txt...";
		dynSimul.saveSensorLogsToFile("LOG.txt");
		cout << "done\n";

		// Show energy balance:
		// -------------------------
		if (SHOW_ENERGY_BALANCE)
		{
			mrpt::gui::CDisplayWindowPlots winE("Energy over time", 600, 300);

			winE.hold_on();
			winE.plot(E_tot, "k3");
			winE.plot(E_kin, "r1");
			winE.plot(E_pot, "b1");

			winE.axis_fit();

			mrpt::gui::CDisplayWindowPlots winE2(
				"Total energy variation", 600, 300);

			vector<double> E_tot_variation = E_tot;
			E_tot_variation += (-E_tot[0]);	 // mrpt::math
			winE2.plot(E_tot_variation, "k2");

			double evMin, evMax;
			mrpt::math::minimum_maximum(E_tot_variation, evMin, evMax);

			winE2.axis(0, E_tot.size(), evMin * 1.1, evMax * 1.1);
			winE2.waitForKey();
		}

		return 0;  // program ended OK.
	}
	catch (exception& e)
	{
		cerr << e.what() << endl;
		return 1;
	}
}
