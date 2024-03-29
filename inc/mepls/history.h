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

#ifndef __HISTORY_H_
#define __HISTORY_H_

#include <mepls/event.h>

#include <deal.II/base/symmetric_tensor.h>

namespace mepls
{

namespace system
{
template<int dim>
class System;
}

/*! @namespace mepls::history 
 * @brief This namespace contains classes and structs used to record the occurrence of
 * events in the system, and ease their write into output files. */
namespace history
{

/*! @class mepls::history::PlasticRow
 * @brief This structs is used by @ref History<dim>::add_row to reformat
* objects of the class @ref event::Plastic as rows of an output dataset.
* The member @ref index is set by @ref History when the event is added. */
struct PlasticRow
{
	unsigned int index;
	unsigned int element;
	double dtime;
	float slip_threshold;
	float dplastic_strain;
	float eigenstrain_00;
	float eigenstrain_11;
	float eigenstrain_01;
	float acting_stress_00;
	float acting_stress_11;
	float acting_stress_01;
	unsigned int activation_protocol;
};

/*! @class mepls::history::DrivingRow
 * @brief This structs is used by @ref History<dim>::add_row to reformat
* objects of the class @ref event::Driving as rows of an output dataset.
* The member @ref index is set by @ref History when the event is added. */
struct DrivingRow
{
	unsigned int index;
	double dext_stress;
	double dpressure;
	double dtotal_strain;
	double dload;
	double dtime;
	unsigned int activation_protocol;
};

/*! @class mepls::history::RenewSlipRow
 * @brief This structs is used by @ref History<dim>::add_row to reformat
* objects of the class @ref event::RenewSlip as rows of an output dataset.
* The member @ref index is set by @ref History when the event is added. */
struct RenewSlipRow
{
	unsigned int index;
	unsigned int element;
	float threshold;
	float slip_angle;
};


/*! @class mepls::history::MacroSummaryRow
 * @brief This structs is used by @ref History<dim>::add_macro to store the system's
 * macroscale propreties. It contains the spatial average and standard deviation of stress tensor
 * components, the von Mises stress, the von Mises plastic strain, the elastic energy, the
 * configurational energy, the time, the total strain, and the external stress
 * */
struct MacroSummaryRow
{
	/*! Struct to store system-scale properties */
	double av_stress_00 = 0.;
	double av_stress_11 = 0.;
	double av_stress_01 = 0.;
	double std_stress_00 = 0.;
	double std_stress_11 = 0.;
	double std_stress_01 = 0.;
	double av_vm_plastic_strain = 0.;
	double std_vm_plastic_strain = 0.;
	double av_vm_stress = 0.;
	double std_vm_stress = 0.;
	double av_energy_el = 0.;
	double std_energy_el = 0.;
	double av_energy_conf = 0.;
	double std_energy_conf = 0.;
	double time = 0.;
	double total_strain = 0.;
	double ext_stress = 0.;
	unsigned int index = 0;
};


/*! class mepls::history::History
 * @brief This class registers the evolution of the system's macroscale properties and the occurrence
 *  of events. The events are recorded in different datasets, depending on whether the event is a
 *  platic or driving one.
 *
 * <hr>
 *
 * Let see an example of how the event datasets are created from a dynamics consisting
 *  in the iterative application of @ref dynamics::extremal_dynamics_step followed by @ref
 *  dynamics::relaxation :
 * - Index=1, driving event: by extremal dynamics, a rise of the external stress
 * - Index=2, plastic event: by relaxation, a mesoscale element undergoes a single slip event after
 *  the previous stress rise
 * - Index=3, driving event: external stress drop due to plastic deformation (under
 * strain-controlled conditions)
 * - Index=4, renew event: modify the structural properties of the deformed element
 * - Index=5, plastic event + plastic event + plastic event +...: by relaxation, several
 * slip eventes due to the stress variations induced by the previous slip event. These
 * events share the same index=5.
 * - Index=6, driving event: external stress drop due to plastic deformation (under
 * strain-controlled conditions)
 * - Index=7, renew event +  renew event + renew event + ...: modify the structural
 * properties of the deformed elements. These events share the same index=7
 * - Index=8, driving event: by extremal dynamics, a rise of the external stress
 * - ...
 *
 * */
template<int dim>
class History
{
  public:

	History(const std::string &inputname_ = "history")
		:
		index_(0),
		closed_(false),
		name_(inputname_)
	{
		/*! Constructor. */
	}

	void close()
	{
		/*! If called, no more data will be recored in the history. */

		closed_ = true;
	}

	void add_row(const event::Driving<dim> &event)
	{
		/*! Dump an object of the class @ref event::Driving into a simple struct
		 * representing a row of an output dataset. */

		DrivingRow row;
		row.dext_stress = event.dext_stress;
		row.dpressure = event.dpressure;
		row.dtotal_strain = event.dtotal_strain;
		row.dload = event.dload;
		row.dtime = event.dtime;
		row.activation_protocol = event.activation_protocol;
		row.index = index_;
		driving.push_back(row);
	}

	void add_row(const event::Plastic<dim> &event)
	{
		/*! Dump an object of the class @ref event::Plastic into a simple struct
		 * representing a row of an output dataset. */

		PlasticRow row;
		row.element = event.element;
		row.dtime = event.dtime;
		row.slip_threshold = event.slip_threshold;
		row.dplastic_strain = event.dplastic_strain;
		row.eigenstrain_00 = float(event.eigenstrain[0][0]);
		row.eigenstrain_11 = float(event.eigenstrain[1][1]);
		row.eigenstrain_01 = float(event.eigenstrain[0][1]);
		row.acting_stress_00 = float(event.acting_stress[0][0]);
		row.acting_stress_11 = float(event.acting_stress[1][1]);
		row.acting_stress_01 = float(event.acting_stress[0][1]);
		row.activation_protocol = event.activation_protocol;
		row.index = index_;
		plastic.push_back(row);
	}

	void add_row(const event::RenewSlip<dim> &event)
	{
		/*! Dump an object of the class @ref event::RenewSlip into a simple struct
		 * representing a row of an output dataset. */

		RenewSlipRow row;
		row.element = event.element_number;
		row.threshold = float(event.threshold);
		row.slip_angle = float(event.slip_angle);
		row.index = index_;
		renew.push_back(row);
	}

	template<typename EventType>
	void add(const EventType &event)
	{
		/*! This function is a (template) wrapper around the different
		 * overloaded functions @ref History.add_row(). It calls
		 * @ref History.add_row() and updates by one unit @ref History.index. */

		if(closed_)
			return;

		add_row(event);
		++index_;
	}

	template<typename EventType>
	void add(const std::vector<EventType> &event_vector)
	{
		/*! The same as @ref History.add(const EventType &event) but
		 * takes a vector of events instead. Those events are added with
		 * the same @ref History.index. */

		if(closed_)
			return;

		for(auto &event : event_vector)
			add_row(event);
		++index_;
	}

	template<typename EventType>
	void add(const std::vector<EventType *> &event_vector)
	{
		/*! The same as @ref History.add(const std::vector<EventType>
		 * &event_vector) but takes a vector of pointers to events. Those events
		 * are added with the same @ref History.index. */

		if(closed_)
			return;

		for(auto &event : event_vector)
			add_row(*event);
		++index_;
	}

	void clear()
	{
		/*! Clear the history. */

		driving.clear();
		plastic.clear();
		renew.clear();
		macro_evolution.clear();
		index_ = 0;
	}


	virtual void add_macro(const mepls::system::System<dim> &system)
	{
		/*! Record the macroscale properties of the input system. By default, it computes
		 * the spatial average and standard deviation of stress tensor components, the von Mises
		 * stress, the von Mises plastic strain, the elastic energy, the configurational
		 * energy, the time, the total strain, and the external stress. */

		if(closed_)
			return;

		MacroSummaryRow data;

		add_macro_default(system, data);

		macro_evolution.push_back(data);
	}


	std::string name() const
	{
		/*! Get the name of the history object. */

		return name_;
	}

	unsigned int index() const
	{
		/*! Get the current history index. */

		return index_;
	}

	std::vector<DrivingRow> driving;
	/*!< Vector with the information of the added driving events.  */

	std::vector<PlasticRow> plastic;
	/*!< Vector with the information of the added plastic events. */

	std::vector<RenewSlipRow> renew;
	/*!< Vector with the information of the microstructural properties renewal events */

	std::vector<MacroSummaryRow> macro_evolution;
	/*!< Vector to store the macroscale properties. */

  protected:


	void add_macro_default(const mepls::system::System<dim> &system, MacroSummaryRow & data)
	{
		/*! Record the default macroscale properties of the input system on the input data struct. */

		const auto &macrostate = system.macrostate;

		double sum_stress_00 = 0.;
		double sum_stress_11 = 0.;
		double sum_stress_01 = 0.;
		double sum_vm_plastic_strain = 0.;
		double sum_vm_stress = 0.;
		double sum_energy_el = 0.;
		double sum_energy_conf = 0.;

		double sum2_stress_00 = 0.;
		double sum2_stress_11 = 0.;
		double sum2_stress_01 = 0.;
		double sum2_vm_plastic_strain = 0.;
		double sum2_vm_stress = 0.;
		double sum2_energy_el = 0.;
		double sum2_energy_conf = 0.;

		for(auto &element : system)
		{
			auto &stress = element->stress();

			sum_stress_00 += stress[0][0];
			sum2_stress_00 += stress[0][0] * stress[0][0];

			sum_stress_11 += stress[1][1];
			sum2_stress_11 += stress[1][1] * stress[1][1];

			sum_stress_01 += stress[0][1];
			sum2_stress_01 += stress[0][1] * stress[0][1];

			double vm_plastic_strain = element->integrated_vm_eigenstrain();
			sum_vm_plastic_strain += vm_plastic_strain;
			sum2_vm_plastic_strain += vm_plastic_strain * vm_plastic_strain;

			double vm_stress = mepls::utils::get_von_mises_equivalent_stress(element->stress());
			sum_vm_stress += vm_stress;
			sum2_vm_stress += vm_stress * vm_stress;

			double energy_el = element->energy_el();
			sum_energy_el += energy_el;
			sum2_energy_el += energy_el * energy_el;

			double energy_conf = element->energy_conf();
			sum_energy_conf += energy_conf;
			sum2_energy_conf += energy_conf * energy_conf;
		}

		double N = double(system.size());

		data.av_stress_00 = sum_stress_00 / N;
		data.av_stress_11 = sum_stress_11 / N;
		data.av_stress_01 = sum_stress_01 / N;
		data.av_vm_stress = sum_vm_stress / N;
		data.av_vm_plastic_strain = sum_vm_plastic_strain / N;
		data.av_energy_el = sum_energy_el / N;
		data.av_energy_conf = sum_energy_conf / N;

		data.std_stress_00 = std::sqrt(sum2_stress_00 / N - data.av_stress_00 * data.av_stress_00);

		data.std_stress_11 = std::sqrt(sum2_stress_11 / N - data.av_stress_11 * data.av_stress_11);

		data.std_stress_01 = std::sqrt(sum2_stress_01 / N - data.av_stress_01 * data.av_stress_01);

		data.std_vm_plastic_strain = std::sqrt(
			sum2_vm_plastic_strain / N - data.av_vm_plastic_strain * data.av_vm_plastic_strain);

		data.std_vm_stress = std::sqrt(sum2_vm_stress / N - data.av_vm_stress * data.av_vm_stress);

		data.std_energy_el = std::sqrt(sum2_energy_el / N - data.av_energy_el * data.av_energy_el);

		data.std_energy_conf = std::sqrt(
			sum2_energy_conf / N - data.av_energy_conf * data.av_energy_conf);

		data.time = macrostate["time"];
		data.total_strain = macrostate["total_strain"];
		data.ext_stress = macrostate["ext_stress"];

		data.index = index_;
	}

	unsigned int index_;
	/*!< This member keeps track of the index to be assigned to the next
	 * registered event. It acts as a global identifier of each
	 * event, independent of its type. Hence, it corresponds to the order
	 * in which each change in the system occurs. The index is updated
	 * one unit in @ref History.add() each time an event is added to
	 * @ref History. Sometimes, it can have the same value for different
	 * events if those events are assumed to occur simultaneously.
	 * For example, if different thresholds within the same element are
	 * simultaneously renewed after plastic deformation occurs in that element.
	 *
	 * @note in the output datasets, it corresponds to the column "index".
	 * It acts as a surrgate key in a relational database that allows
	 * reconstructing the simulation history by merging the datasets. */

	bool closed_;
	/*!< If true, no further events are added. */

	std::string name_;
	/*!< The name of the history object. */
};


} // namespace history


} // namespace mepls

#endif //__HISTORY_H_
