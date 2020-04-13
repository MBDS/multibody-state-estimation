
// Example of MBS model: Double Quad
// -------------------------------------
#include <mbse/sparsembs.h>
#include <mbse/model-examples.h>

#include <mrpt/opengl.h>
#include <mrpt/gui.h>
#include <mrpt/math/ops_vectors.h>
#include <thread>  // for sleep()

using namespace std;
using namespace mbse;
using namespace mrpt;
using namespace mrpt::poses;
using namespace mrpt::math;

#define SIMUL_REALTIME 1  // 1: real time, 0: fixed simulation step

double REALTIME_FACTOR = 0.15;

const double GUI_DESIRED_FPS = 75;  // Hz

const bool SHOW_ENERGY_BALANCE = true;

void my_callback(const TSimulationStateRef simul_state) {}

int main(int argc, char** argv)
{
	try
	{
		CBody::TRenderParams dynamic_rp;

		CModelDefinition model;

		// Define a double pendulum mechanism:
		// ------------------------------------------------
		buildFourBarsMBS(model);
		// buildFourBarsMBS_JavierCuadrado(model);
		// buildSliderCrankMBS(model);
		// buildFollowerMBS(model);
		// buildLongStringMBS( 15, model);
		// sparsembs::buildTwoSliderBlocks(model);

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
		// CDynamicSimulator_Indep_dense
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

		mrpt::opengl::TFontParams fp;
		fp.vfont_scale = 15;  // pixels
		fp.draw_shadow = true;

		while (win3D.isOpen())
		{
			double t_new = tictac.Tac();

#if SIMUL_REALTIME
			if (t_new >
				t_old + dynSimul.params.time_step)  // Just in case the computer
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
					0 /* txt ID */, fp);

				const double simul_t =
					sparsembs::timelog.getMeanTime("mbs.run_complete_timestep");
				const double simul_Hz = simul_t > 0 ? 1.0 / simul_t : 0;
				win3D.addTextMessage(
					10, 30,
					mrpt::format(
						"Simul: %.02fHz [%.02fms]", simul_Hz, simul_t * 1e3),
					1 /* txt ID */, fp);

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
			E_tot_variation += (-E_tot[0]);  // mrpt::math
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
