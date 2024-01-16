/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2024 University of Almeria                           |
  | Copyright (C) 2021 University of Salento                                |
  | See README for list of authors and papers                               |
  | Distributed under 3-clause BSD license                                  |
  |  See: <https://opensource.org/licenses/BSD-3-Clause>                    |
  +-------------------------------------------------------------------------+ */

// Loads a mechanism and render it and, optionally, a trajectory.

#include <mrpt/core/lock_helper.h>
#include <mrpt/core/round.h>
#include <mbse/mbse.h>
#include <mrpt/opengl.h>
#include <mrpt/gui/CDisplayWindowGUI.h>
//#include <mrpt/math/ops_vectors.h>
#include <thread>  // for sleep()
#include <mrpt/system/filesystem.h>

#include <mrpt/3rdparty/tclap/CmdLine.h>

using namespace std;
using namespace mbse;
using namespace mrpt;
using namespace mrpt::poses;
using namespace mrpt::math;

TCLAP::CmdLine cmd("mbse-viewer", ' ');

TCLAP::ValueArg<std::string> arg_mechanism(
	"", "mechanism", "Mechanism model YAML file", true, "mechanism.yaml",
	"YAML model definition", cmd);

TCLAP::ValueArg<std::string> arg_q(
	"q", "q", "Trajectory 'q' (positions) TXT file", false, "q.txt", "q.txt",
	cmd);

TCLAP::ValueArg<std::string> arg_dq(
	"", "dq", "Trajectory 'dq' (velocities) TXT file", false, "dq.txt",
	"dq.txt", cmd);

TCLAP::ValueArg<std::string> arg_ddq(
	"", "ddq", "Trajectory 'ddq' (accelerations) TXT file", false, "ddq.txt",
	"ddq.txt", cmd);

static void runViewer()
{
	// Load mechanism model:
	const auto yamlData =
		mrpt::containers::yaml::FromFile(arg_mechanism.getValue());
	const ModelDefinition model = ModelDefinition::FromYAML(yamlData);

	// "Compile" the problem:
	AssembledRigidModel::Ptr aMBS = model.assembleRigidMBS();

	std::cout << "Problem coordinates:\n";
	aMBS->printCoordinates(std::cout);
	std::cout << std::endl;

	if (0)
	{
		std::cout << "Problem constraints:\n";
		aMBS->printConstraints(std::cout);
		std::cout << std::endl;
	}

	mrpt::math::CMatrixDouble q, dq, ddq;
	if (arg_q.isSet()) q.loadFromTextFile(arg_q.getValue());

	const size_t nTimeSteps = q.rows();
	if (nTimeSteps > 0 && arg_dq.isSet())
	{
		dq.loadFromTextFile(arg_dq.getValue());
		ASSERT_EQUAL_(q.rows(), dq.rows());
		ASSERT_EQUAL_(q.cols(), dq.cols());
	}
	if (nTimeSteps > 0 && arg_ddq.isSet())
	{
		ddq.loadFromTextFile(arg_ddq.getValue());
		ASSERT_EQUAL_(q.rows(), ddq.rows());
		ASSERT_EQUAL_(q.cols(), ddq.cols());
	}
	// Remove first time columns:
	std::vector<double> timestamps;
	if (nTimeSteps > 0)
	{
		q.extractColumn(0, timestamps);

		q.removeColumns({0});
		if (dq.cols() > 0) dq.removeColumns({0});
		if (ddq.cols() > 0) ddq.removeColumns({0});
	}

	// Prepare 3D scene:
	// -----------------------------------------------
	Body::TRenderParams dynamic_rp;
	auto gl_MBS = mrpt::opengl::CSetOfObjects::Create();
	aMBS->getAs3DRepresentation(gl_MBS, dynamic_rp);

	nanogui::init();

	mrpt::gui::CDisplayWindowGUI win(
		mrpt::format(
			"mbse-viewer [%s]",
			mrpt::system::extractFileName(arg_mechanism.getValue()).c_str()),
		1280, 768);

	auto scene = mrpt::opengl::COpenGLScene::Create();

	{
		auto lck = mrpt::lockHelper(win.background_scene_mtx);
		auto& cam = win.camera();
		cam.setAzimuthDegrees(-90);
		cam.setElevationDegrees(90);
		cam.setZoomDistance(15);

		// gl_MBS->setPose(mrpt::poses::CPose3D(0, 0, 0, DEG2RAD(180),
		// DEG2RAD(0), DEG2RAD(90)));
		scene->insert(mrpt::opengl::stock_objects::CornerXYZSimple(0.1));
		scene->insert(gl_MBS);

		win.background_scene = scene;
	}

	// subwindow:
	auto w = win.createManagedSubWindow("Mechanism");
	w->setPosition({5, 80});
	w->setLayout(new nanogui::BoxLayout(
		nanogui::Orientation::Vertical, nanogui::Alignment::Fill));
	w->setFixedWidth(260);

	// time slider:
	auto lbTime = w->add<nanogui::Label>("Time: [no trajectory loaded]");
	auto slider = w->add<nanogui::Slider>();
	slider->setRange({0, 1});
	if (nTimeSteps == 0) slider->setEnabled(false);

	// Tabs:
	auto tab = w->add<nanogui::TabWidget>();

	auto lbPhi = w->add<nanogui::Label>("|Phi|=0");
	auto btnPosProblem = w->add<nanogui::Button>("Position problem");

	std::vector<nanogui::Widget*> tabs = {
		tab->createTab("q"), tab->createTab("dq"), tab->createTab("ddq")};

	tab->setActiveTab(0);

	for (auto t : tabs)
		t->setLayout(new nanogui::BoxLayout(
			nanogui::Orientation::Vertical, nanogui::Alignment::Fill, 1, 1));

	const int pnWidth = 240, pnHeight = 270;

	std::vector<nanogui::VScrollPanel*> vscrolls;
	for (auto t : tabs) vscrolls.emplace_back(t->add<nanogui::VScrollPanel>());

	// vscroll should only have *ONE* child.
	// this is what `wrapper` is for
	std::vector<nanogui::Widget*> wrappers;

	for (auto vs : vscrolls)
	{
		vs->setFixedSize({pnWidth, pnHeight});
		auto wr = vs->add<nanogui::Widget>();
		wr->setLayout(new nanogui::GridLayout(
			nanogui::Orientation::Horizontal, 2 /*columns */,
			nanogui::Alignment::Minimum, 2, 2));

		wr->setFixedSize({pnWidth - 7, pnHeight});
		wrappers.emplace_back(wr);
	}

	// q, dq, ddq
	std::vector<nanogui::TextBox*> tbs[3];
	std::vector<nanogui::CheckBox*> cbs[3];

	for (int pnIdx = 0; pnIdx < 3; pnIdx++)
	{
		auto pn = wrappers.at(pnIdx);

		const char* varName;
		switch (pnIdx)
		{
			case 0:
				varName = "q";
				break;
			case 1:
				varName = "dq";
				break;
			case 2:
				varName = "ddq";
				break;
		};

		for (int i = 0; i < aMBS->q_.size(); i++)
		{
			auto cb =
				pn->add<nanogui::CheckBox>(mrpt::format("%s[%i]", varName, i));
			cbs[pnIdx].push_back(cb);

			auto tb = pn->add<nanogui::TextBox>();
			tb->setEditable(true);
			tbs[pnIdx].push_back(tb);
		}
	}

	const auto lambdaUpdateQValues = [&]() {
		for (int pnIdx = 0; pnIdx < 3; pnIdx++)
		{
			Eigen::VectorXd v;
			switch (pnIdx)
			{
				case 0:
					v = aMBS->q_;
					break;
				case 1:
					v = aMBS->dotq_;
					break;
				case 2:
					v = aMBS->ddotq_;
					break;
			};

			for (int i = 0; i < v.size(); i++)
				tbs[pnIdx][i]->setValue(mrpt::format("%.3f", v[i]));
		}

		// Indicators:
		lbPhi->setCaption(mrpt::format("|Phi|=%g", aMBS->Phi_.norm()));

		// refresh 3D objects:
		{
			auto lck = mrpt::lockHelper(win.background_scene_mtx);
			aMBS->update3DRepresentation(dynamic_rp);
		}
	};

	lambdaUpdateQValues();

	if (nTimeSteps > 0)
	{
		slider->setRange({0, nTimeSteps - 1});
		slider->setCallback([&](float pos) {
			const auto timIdx = mrpt::round(pos);
			lbTime->setCaption(
				mrpt::format("Time: %10.05f s", timestamps.at(timIdx)));

			q.extractRow(timIdx, aMBS->q_);
			if (dq.rows() > 0) dq.extractRow(timIdx, aMBS->dotq_);
			if (ddq.rows() > 0) ddq.extractRow(timIdx, aMBS->ddotq_);

			aMBS->update_numeric_Phi_and_Jacobians();

			lambdaUpdateQValues();
		});

		slider->setValue(0);
		slider->callback()(slider->value());
	}

	btnPosProblem->setCallback([&]() {
		std::vector<size_t> idxs;
		for (size_t i = 0; i < cbs[0].size(); i++)
		{
			if (cbs[0][i]->checked())
			{
				idxs.push_back(i);
				const double v = std::stod(tbs[0][i]->value());
				aMBS->q_[i] = v;
			}
		}

		AssembledRigidModel::ComputeDependentParams cdp;
		AssembledRigidModel::ComputeDependentResults cdr;

		aMBS->computeDependentPosVelAcc(idxs, true, false, cdp, cdr);

		lambdaUpdateQValues();
	});

	//
	win.performLayout();

	win.drawAll();
	win.setVisible(true);
	nanogui::mainloop();

	nanogui::shutdown();
}

int main(int argc, char** argv)
{
	try
	{
		// Parse arguments:
		if (!cmd.parse(argc, argv))
			throw std::runtime_error("");  // should exit.

		runViewer();

		return 0;  // program ended OK.
	}
	catch (exception& e)
	{
		cerr << e.what() << endl;
		return 1;
	}
}
