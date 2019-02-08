#pragma once

#include "sparsembs-common.h"
#include "sparsembs-utils.h"
#include "CAssembledModelRigid.h"

namespace sparsembs
{
using namespace Eigen;
using namespace mrpt::math;

/** Base of all types of virtual sensors */
class CVirtualSensor
{
   public:
	/** Simulates one sensor reading from the given system state, returning the
	 * predicted value. */
	virtual double simulate_reading(
		const CAssembledRigidModel& mb_state) const = 0;

	/** Returns the log-likelihhod of the given read value for the current
	 * mechanism state */
	double evaluate_log_likelihood(
		const double sensor_reading, const CAssembledRigidModel& mb_state) const
	{
		const double sensor_prediction = this->simulate_reading(mb_state);
		return -0.5 *
			   mrpt::utils::square(
				   (sensor_reading - sensor_prediction) / sensor_noise_std);
	}

	/** One standard deviation (1sigma) of the sensor Gaussian noise model
	 * (units are sensor-specific) */
	double sensor_noise_std;

	CVirtualSensor() : sensor_noise_std(1) {}

	virtual ~CVirtualSensor() {}
};

typedef stlplus::smart_ptr<CVirtualSensor> CVirtualSensorPtr;

/** A Gyroscope sensor */
class CVirtualSensor_Gyro : public CVirtualSensor
{
   public:
	/** Simulates one sensor reading from the given system state, returning the
	 * predicted value. */
	virtual double simulate_reading(const CAssembledRigidModel& mb_state) const
		MRPT_OVERRIDE;

	CVirtualSensor_Gyro(const size_t body_idx) : m_body_idx(body_idx) {}

   protected:
	size_t m_body_idx;
};

}  // namespace sparsembs
