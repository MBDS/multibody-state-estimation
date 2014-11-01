#pragma once

#include "sparsembs-common.h"
#include "sparsembs-utils.h"
#include "CBody.h"
#include "constraints.h"
#include "CModelDefinition.h"
#include "CAssembledModelRigid.h"
#include "dynamic-simulators.h"
#include "virtual-sensors.h"

#include <mrpt/bayes/CParticleFilter.h>
#include <mrpt/bayes/CParticleFilterCapable.h>
#include <mrpt/bayes/CParticleFilterData.h>
#include <mrpt/random.h>

namespace sparsembs
{
	using namespace Eigen;
	using namespace mrpt::math;

	/** A "particle" holding a MultiBody assembly + its state */
	template <class DYN_SIMUL>
	struct TMBState_Particle
	{
		TMBState_Particle(const TSymbolicAssembledModel &sym_model) :
			sym_assembly_model(sym_model),
			num_model_ptr(new CAssembledRigidModel(sym_assembly_model)),
			num_model(*num_model_ptr.pointer()),
			dyn_simul(new DYN_SIMUL(num_model_ptr))
		{
			dyn_simul->prepare();
		}

		TSymbolicAssembledModel  sym_assembly_model; //!< Symbolic assembly instructions
		CAssembledRigidModelPtr  num_model_ptr;
		CAssembledRigidModel &   num_model;
		CDynamicSimulatorIndepBasePtr dyn_simul;

		/** copy ctor */
		TMBState_Particle(const TMBState_Particle &o) :
			sym_assembly_model(o.sym_assembly_model),
			num_model_ptr(new CAssembledRigidModel(sym_assembly_model)),
			num_model(*num_model_ptr.pointer()),
			dyn_simul(new DYN_SIMUL(num_model_ptr))
		{
			num_model.copyStateFrom(o.num_model);
			dyn_simul->prepare();
		}

		/** copy operator: doesn't copy the entire objects, just the state! */
		TMBState_Particle & operator =(const TMBState_Particle &o)
		{
			num_model.copyStateFrom(o.num_model);
			return *this;
		}
	};

	// select the kind of solver:
	typedef CDynamicSimulator_Indep_dense  MBPF_SIMULATOR_TYPE;

	/** A particle-based representation of a probability density function (PDF) over the state of a mechanical system.
	  */
	class CMultiBodyParticleFilter
		:
		public mrpt::bayes::CParticleFilterData<TMBState_Particle<MBPF_SIMULATOR_TYPE> >,
		public mrpt::bayes::CParticleFilterDataImpl<CMultiBodyParticleFilter,mrpt::bayes::CParticleFilterData<TMBState_Particle<MBPF_SIMULATOR_TYPE> >::CParticleList>
	{
	public:
		typedef TMBState_Particle<MBPF_SIMULATOR_TYPE> particle_t;


		/** Initializes a set of M particles for the given multibody system */
		CMultiBodyParticleFilter(const size_t M, const CModelDefinition &mbs);

		/** Dtor */
		~CMultiBodyParticleFilter();


		struct TOutputInfo
		{
			bool resampling_done; //!< =true if resampling was required
			double ESS;

			TOutputInfo() : resampling_done(false),ESS(1)
			{ }
		};

		/** Runs one step of the PF (SIR) algorithm.
		  * Simulation runs for "t_increment", but several steps are runned if that value is greater than "max_t_step".
		*/
		void run_PF_step(
			const double t_ini,
			const double t_end,
			const double max_t_step,
			const std::vector<CVirtualSensorPtr> & sensor_descriptions,
			const std::vector<double>          & sensor_readings,
			TOutputInfo &out_info
			);

		void getAs3DRepresentation( mrpt::opengl::CSetOfObjectsPtr & 	outObj,const CBody::TRenderParams & rp ) const;
		void update3DRepresentation(const CBody::TRenderParams & rp) const;



		struct TTransitionModelOptions
		{
			double acc_xy_noise_std; //!< 1 sigma of the additive Gaussian noise for accelerations in X,Y.

			TTransitionModelOptions();
		};

		TTransitionModelOptions model_options;
		mrpt::bayes::CParticleFilter::TParticleFilterOptions PF_options; //!< Parameters for the PF algorithm.

		mrpt::random::CRandomGenerator  random_generator;

	}; // end class CMultiBodyParticleFilter

}

