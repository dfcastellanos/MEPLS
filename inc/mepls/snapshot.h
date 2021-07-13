// -----------------------------------------------------------------------
//
// Copyright (C) 2020  - David Fernández Castellanos
//
// This file is part of the MEPLS software. You can use it, redistribute
// it, and/or modify it under the terms of the Creative Commons Attribution
// 4.0 International Public License. The full text of the license can be
// found in the file LICENSE at the top level of the MEPLS distribution.
//
// -----------------------------------------------------------------------

#ifndef RUN_SIM_SNAPSHOTS_H
#define RUN_SIM_SNAPSHOTS_H

#include <deal.II/base/symmetric_tensor.h>


namespace mepls
{

/*! This namespace contains classes for taking snapshots of the state of the 
 * elements composing the system. A snapshot consists of a copy of one of the 
 * mesoscale fields over all the elements of the system into a vector.
 * Moreover, the snapshots store metadata such as the name of the recorded 
 * field, the value of a monitor magnitude at which they were taken (e.g., time)
 * etc. */
namespace snapshot
{

/*! Snapshot of @ref element::Element.stress. */
template<int dim>
class Stress
{

  public:

	/*! This struct is used for output purposes only. It converts a complex
	 * tensor object into a simple struct of scalar values. In this way, it can
	 * be easily written into output datasets. */
	struct DataRow
	{
		float xx = 0.;
		float yy = 0.;
		float xy = 0.;
	};

	Stress(
		const system::System<dim> &system,
		std::string monitor_mag_,
		double desired_target_,
		double recorded_target_)
		:
		recorded_mag("stress"),
		monitor_name(monitor_mag_),
		desired_target(desired_target_),
		recorded_target(recorded_target_),
		output_index(system.history->index())
	{
		/*! Take and store the snapshot.
		 *
		 * @param system system from where the snapshot is taken
		 * @param monitor_mag_ name of the magnitude used to check whether the
		 * snapshot should be taken or not
		 * @param desired_target_ value of the @ref monitor_name at which we
		 * desired to take the snapshot
		 * @param recorded_target_ value of the @ref monitor_name at which the
		 * snapshot is actually taken */

		for(auto &element : system)
		{ M_Assert(element->number() == data.size(), "Element data not written in the right order");

			DataRow row;
			row.xx = element->stress()[0][0];
			row.yy = element->stress()[1][1];
			row.xy = element->stress()[0][1];
			data.push_back(row);
		}
	};

	std::vector<DataRow> data;
	/*!< Container to store the copied field. */

	std::string recorded_mag;
	/*!< Name of the field storaged in @ref data. */

	std::string monitor_name;
	/*!< Name of the magnitude used to check whether the snapshot should be taken
	 * or not. */

	double desired_target;
	/*!< Value of the @ref monitor_name at which we desired to take the
	 * snapshot. */

	double recorded_target;
	/*!< Value of the @ref monitor_name at which the snapshot is actually
	 * taken. */

	unsigned int output_index;
	/*!< Global event index from the history (see @ref mepls::history::History<dm>)
	 * at which the snapshot is taken. */

	typename std::vector<DataRow>::iterator begin()
	{
		return data.begin();
	};

	typename std::vector<DataRow>::iterator end()
	{
		return data.end();
	};
};


/*! Snapshot of @ref element::Element.def_grad. */
template<int dim>
class DefGrad
{
  public:

	/*! This struct is used for output purposes only. It converts a complex
	 * tensor object into a simple struct of scalar values. In this way, it can
	 * be easily written into output datasets. */
	struct DataRow
	{
		float xx = 0.;
		float yy = 0.;
		float xy = 0.;
		float yx = 0.;
	};

	DefGrad(
		const system::System<dim> &system,
		std::string monitor_mag_,
		double desired_target_,
		double recorded_target_)
		:
		recorded_mag("def_grad"),
		monitor_name(monitor_mag_),
		desired_target(desired_target_),
		recorded_target(recorded_target_),
		output_index(system.history->index())
	{
		/*! Take and store the snapshot.
		 *
		 * @param system system from where the snapshot is taken
		 * @param monitor_mag_ name of the magnitude used to check whether the
		 * snapshot should be taken or not
		 * @param desired_target_ value of the @ref monitor_name at which we
		 * desired to take the snapshot
		 * @param recorded_target_ value of the @ref monitor_name at which the
		 * snapshot is actually taken */

		for(auto &element : system)
		{ M_Assert(element->number() == data.size(), "Element data not written in the right order");

			DataRow row;
			row.xx = element->def_grad()[0][0];
			row.yy = element->def_grad()[1][1];
			row.xy = element->def_grad()[0][1];
			row.yx = element->def_grad()[1][0];
			data.push_back(row);
		}
	};

	std::vector<DataRow> data;
	/*!< Container to store the copied field. */

	std::string recorded_mag;
	/*!< Name of the field storaged in @ref data. */

	std::string monitor_name;
	/*!< Name of the magnitude used to check whether the snapshot should be taken
	 * or not. */

	double desired_target;
	/*!< Value of the @ref monitor_name at which we desire to take the
	 * snapshot. */

	double recorded_target;
	/*!< Value of the @ref monitor_name at which the snapshot is actually
	 * taken. */

	unsigned int output_index;
	/*!< Global event index from the history (see @ref mepls::history::History<dm>)
	 * at which the snapshot is taken. */

	typename std::vector<DataRow>::iterator begin()
	{
		return data.begin();
	};

	typename std::vector<DataRow>::iterator end()
	{
		return data.end();
	};
};


/*! Snapshot of all the slip thresholds (see @ref mepls::slip::Slip<dim>::threshold)
 * in the system. */
template<int dim>
class Threshold
{

  public:

	/*! This struct is used for output purposes only. It converts a complex
	 * @ref mepls::slip::Slip<dim> object into a simple struct of scalar values. In this
	 * way, it can be easily written into output datasets. */
	struct DataRow
	{
		unsigned int element = 0;
		float angle = 0.;
		float threshold = 0.;
	};

	Threshold(
		const system::System<dim> &system,
		std::string monitor_mag_,
		double desired_target_,
		double recorded_target_)
		:
		recorded_mag("threshold"),
		monitor_name(monitor_mag_),
		desired_target(desired_target_),
		recorded_target(recorded_target_),
		output_index(system.history->index())
	{
		/*! Take and store the snapshot.
		 *
		 * @param system system from where the snapshot is taken
		 * @param monitor_mag_ name of the magnitude used to check whether the
		 * snapshot should be taken or not
		 * @param desired_target_ value of the @ref monitor_name at which we
		 * desired to take the snapshot
		 * @param recorded_target_ value of the @ref monitor_name at which the
		 * snapshot is actually taken */

		for(auto &element : system)
			for(auto &slip : *element)
			{
				DataRow row;
				row.element = element->number();
				row.angle = slip->angle;
				row.threshold = slip->threshold;
				data.push_back(row);
			}
	};

	std::vector<DataRow> data;
	/*!< Container to store the slip thresholds followed and slip angles.*/

	std::string recorded_mag;
	/*!< Name of the field storaged in @ref data. */

	std::string monitor_name;
	/*!< Name of the magnitude used to check whether the snapshot should be taken
	 * or not. */

	double desired_target;
	/*!< Value of the @ref monitor_name at which we desired to take the
	 * snapshot. */

	double recorded_target;
	/*!< Value of the @ref monitor_name at which the snapshot is actually
	 * taken. */

	unsigned int output_index;
	/*!< Global event index from the history (see @ref mepls::history::History<dm>)
	 * at which the snapshot is taken. */

	typename std::vector<DataRow>::iterator begin()
	{
		return data.begin();
	};

	typename std::vector<DataRow>::iterator end()
	{
		return data.end();
	};
};


/*! This class decides whether a snapshot should be taken or not based on the
 * value of a certain monitor magnitude, which is compared with the values of an
 * internal list. The values of the list are those values at which we desire to
 * take a snapshot. The comparison of the value of type double is made within 
 * a certain sensitivity. */
template<int dim>
class Check
{
  public:
	Check(double min, double max, double interval, double sensitivity)
	{
		/*! Constructure. Create a list of value from `min` to `max`, in intervals of
		 * length `interval`. The sensitivity of is set to `sensitivity`. */

		double v = min;
		while(v < max)
		{
			v += interval;
			target_values.push_back(v);
		}

		this->sensitivity = sensitivity;
		already_recorded.assign(target_values.size(), false);
	}

	bool operator()(double val)
	{
		/*! Return true if a snapshot must be taken, false otherwise. A snapshot
		 * must be taken if the input value matches, within a tolerance given by
		 * @ref sensitivity, one of the values in @ref target_values. When there
		 * is a match true is returned, and true is set in the corresponding index
		 * of @ref already_recorded, so no future snapshots will be taken again
		 * with the same value (this can happen if future input values still
		 * match the same target value with the given @ref sensitivity).
		 *
		 * @note to which magnitude (e.g., plastic strain, time, etc.) the input
		 * value corresponds does not matter to this class. */

		for(unsigned int i = 0; i < target_values.size(); ++i)
			if((std::abs(val - target_values[i]) < sensitivity) and not already_recorded[i])
			{
				already_recorded[i] = true;
				desired_value = target_values[i];
				return true;
			}

		return false;
	}

	double sensitivity;
	/*!< Numerical tolerance to decide whether a double matches one of the double
	 * values in @ref target_values. */

	std::vector<double> target_values;
	/*!< List of values were at which snapshots should be taken. */

	std::vector<bool> already_recorded;
	/*!< Mark as true the target values that have already been matched to avoid
	 * taking snapshot again at the same value. */

	double desired_value;
	/*!< Value at which we desire to take the snapshot (thus it is one in @ref
	 * target_values). */
};

} // namespace snapshot
} // namespace mepls

#endif //RUN_SIM_SNAPSHOTS_H
