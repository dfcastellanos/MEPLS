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

#include <deal.II/base/symmetric_tensor.h>


namespace mepls
{

namespace slip
{
template<int dim>
class Slip;
}

namespace system
{
template<int dim>
class System;
}

/*! This namespace contains Event classes, which are meant to represent
 * the occurrence of different kinds of discrete events within the system.
 * It contains as well as the History class which is used to register the
 * occurrence of the events. */
namespace event
{

/*! This class represents a local plastic deformation, which is the consequence 
 * of a slip system activation event */
template<int dim>
struct Plastic
{
	Plastic(slip::Slip<dim> *slip_)
		:
		slip(slip_),
		dplastic_strain(0.),
		dtime(0.),
		acting_stress(slip_->parent->stress()),
		activation_protocol(0),
		element(slip_->parent->number())
	{
		/*! Constructor.
		 * @warning the input pointer can be invalidated if the slip's parent
		 * element \ref slip::Slip::parent gets its structural properties renewed
		 * by calling to element::Element::renew_structural_properties, and the
		 * implementation of such function by derived element classes involves
		 * removing the existing slips. */
	};

	Plastic()
		:
		dplastic_strain(0.),
		dtime(0.),
		activation_protocol(0)
	{
		/*! Constructor.
		 *
		 * @warning after constructing this object, the member \ref slip remains
		 * uninitialized. */
	};

	slip::Slip<dim> *slip;
	/*!<  Pointer to the Slip object which has been activated. */

	double dplastic_strain;
	/*!< Scalar effective plastic shear strain (see \ref
	 * element::Element.vm_eigenstrain). */

	double dtime;
	/*!< Duration of the event. */

	dealii::SymmetricTensor<2, dim> eigenstrain;
	/*!< Eigenstrain added to the parent element (see \ref base::Slip.parent) of
	 * the activated slip object due to the plastic event. See \ref
	 * element::Element.eigenstrain. */

	dealii::SymmetricTensor<2, dim> acting_stress;
	/*!< Stress tensor acting on the parent element (see \ref base::Slip.parent)
	 * at the moment of the platic event. */

	unsigned int activation_protocol;
	/*!< Dynamics protocol (see \ref dynamics::Protocol) by which the plastic
	 * event has been triggered. */

	unsigned int element;
	/*!< Number of the element in which the plastic deformation takes place. */
};


/*! This class represents changes in the external driving conditions. */
template<int dim>
struct Driving
{
	Driving()
		:
		dext_stress(0.),
		dpressure(0.),
		dtotal_strain(0.),
		dload(0.),
		dtime(0.),
		activation_protocol(0)
	{
		/*! Constructor. Initialize members to 0. */
	};

	double dpressure;

	double dext_stress;
	/*!< Change of the applied external stress. It is computed in \ref
	 * System.update_stress_field(). A change of external stress occurs when the
	 * external load varies or when plastic deformation occurs in
	 * strain-controlled conditions. */

	double dtotal_strain;
	/*!< Change of the total strain (i.e., the sum of plastic strain and
	 * elastic strain). It is computed in \ref System.update_stress_field(). */

	double dload;
	/*!< Change of the external load. */

	double dtime;
	/*!< Time interval over which the change of driving conditions occurs. */

	unsigned int activation_protocol;
	/*!< Dynamics protocol (see \ref dynamics::Protocol) responsible of
	 * varying the external driving conditions.*/
};


/*! This class represents changes in the local microstructural properties. */
template<int dim>
struct RenewSlip
{
	RenewSlip(slip::Slip<dim> *slip)
		:
		element_number(slip->parent->number()),
		threshold(slip->threshold),
		slip_angle(slip->angle)
	{
		/*! Constructor. Initialize the pointers. */
	};

	unsigned int element_number;
	/*!< Number of the parent element of the Slip object. */

	double threshold;
	/*!< Copy of \ref Slip.threshold. */

	double slip_angle;
	/*!< Copy of \ref Slip.angle. */
};

} // namespace event


/*! This namespace contains classes and structs used to record the occurrence of
 * events in the system, and ease their write into output files. */
namespace history
{

/*! This structs is used by \ref History.add_row() to reformat
* objects of the class \ref event::Plastic as rows of an output dataset.
* The member \ref index is set by \ref History when the event is added. */
struct PlasticRow
{
	unsigned int index;
	unsigned int element;
	double dtime;
	float eigenstrain_00;
	float eigenstrain_11;
	float eigenstrain_01;
	float acting_stress_00;
	float acting_stress_11;
	float acting_stress_01;
	unsigned int activation_protocol;
};

/*! This structs is used by \ref History.add_row() to reformat
* objects of the class \ref event::Driving as rows of an output dataset.
* The member \ref index is set by \ref History when the event is added. */
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

/*! This structs is used by \ref History.add_row() to reformat
* objects of the class \ref event::RenewSlip as rows of an output dataset.
* The member \ref index is set by \ref History when the event is added. */
struct RenewSlipRow
{
	unsigned int index;
	unsigned int element;
	float threshold;
	float slip_angle;
};

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
	double av_potential_energy = 0.;
	double std_potential_energy = 0.;
	double time = 0.;
	double total_strain = 0.;
	double ext_stress = 0.;
	unsigned int index = 0;
};


/*! This class registers the occurrence of events. The data contained in
 * the event objects are reformated into structs in a way that makes easier
 * the simulation output. Each event is indexed by a global index simply
 * referred to as “index”, which denotes the number of events that occurred in
 * history object, irrespective of the type of event.
 *
 * <hr>
 *
 * Let see an example that results from a dynamics consisting in the iterative
 * application of protocol 2 called (see \ref dynamics::Protocol::extremal_dynamics
 * and \ref dynamics::extremal_dynamics_step), followed by protocol 4
 * (see \ref dynamics::Protocol::relaxation_step and \ref dynamics::relaxation):
 * - Index=1, driving event: by extremal dynamics, a rise of the external stress
 * - Index=2, plastic event: by relaxation, a mesoscale element yields after the
 * previous stress rise
 * - Index=3, driving event: external stress drop due to plastic deformation under
 * strain-controlled conditions
 * - Index=4, renew event: modify the structural properties of the yielded element
 * - Index=5, plastic event + plastic event + plastic event: by relaxation, 3
 * different elements simultaneously yield due to the stress variations induced
 * by the previous plastic event. These three events share the same index=5.
 * - Index=6, driving event: external stress drop due to plastic deformation under
 * strain-controlled conditions
 * - Index=7, renew event +  renew event + renew event: modify the structural
 * properties of the yielded elements. These three events share the same index=7
 * - Index=8, driving event: by extremal dynamics, a rise of the external stress
 * - ...
 *
 * */
template<int dim>
class History
{
  public:

	History(const std::string &name_ = "history")
		:
		index(0),
		closed(false),
		name(name_)
	{
		/*! Constructor. */
	}

	void close()
	{
		closed = true;
	}

	void add_row(const event::Driving<dim> &event)
	{
		/*! Dump an object of the class \ref event::Driving into a simple struct
		 * representing a row of an output dataset. */

		DrivingRow row;
		row.dext_stress = event.dext_stress;
		row.dpressure = event.dpressure;
		row.dtotal_strain = event.dtotal_strain;
		row.dload = event.dload;
		row.dtime = event.dtime;
		row.activation_protocol = event.activation_protocol;
		row.index = index;
		driving.push_back(row);
	}

	void add_row(const event::Plastic<dim> &event)
	{
		/*! Dump an object of the class \ref event::Plastic into a simple struct
		 * representing a row of an output dataset. */

		PlasticRow row;
		row.element = event.element;
		row.dtime = event.dtime;
		row.eigenstrain_00 = float(event.eigenstrain[0][0]);
		row.eigenstrain_11 = float(event.eigenstrain[1][1]);
		row.eigenstrain_01 = float(event.eigenstrain[0][1]);
		row.acting_stress_00 = float(event.acting_stress[0][0]);
		row.acting_stress_11 = float(event.acting_stress[1][1]);
		row.acting_stress_01 = float(event.acting_stress[0][1]);
		row.activation_protocol = event.activation_protocol;
		row.index = index;
		plastic.push_back(row);
	}

	void add_row(const event::RenewSlip<dim> &event)
	{
		/*! Dump an object of the class \ref event::RenewSlip into a simple struct
		 * representing a row of an output dataset. */

		RenewSlipRow row;
		row.element = event.element_number;
		row.threshold = float(event.threshold);
		row.slip_angle = float(event.slip_angle);
		row.index = index;
		renew.push_back(row);
	}

	template<typename EventType>
	void add(const EventType &event)
	{
		/*! This function is a (template) wrapper around the different
		 * overloaded functions \ref History.add_row(). It calls
		 * \ref History.add_row() and updates by one unit \ref History.index. */

		if(closed)
			return;

		add_row(event);
		++index;
	}

	template<typename EventType>
	void add(const std::vector<EventType> &event_vector)
	{
		/*! The same as \ref History.add(const EventType &event) but
		 * takes a vector of events instead. Those events are added with
		 * the same \ref History.index. */

		if(closed)
			return;

		for(auto &event : event_vector)
			add_row(event);
		++index;
	}

	template<typename EventType>
	void add(const std::vector<EventType *> &event_vector)
	{
		/*! The same as \ref History.add(const std::vector<EventType>
		 * &event_vector) but takes a vector of pointers to events. Those events
		 * are added with the same \ref History.index. */

		if(closed)
			return;

		for(auto &event : event_vector)
			add_row(*event);
		++index;
	}

	void clear()
	{
		driving.clear();
		plastic.clear();
		renew.clear();
		macro_evolution.clear();
		index = 0;
	}


   virtual void add_macro(const mepls::system::System<dim> &system)
   {
      if(closed)
         return;

		MacroSummaryRow data;

		const auto &elements = system.elements;
		const auto &macrostate = system.macrostate;

		double sum_stress_00 = 0.;
		double sum_stress_11 = 0.;
		double sum_stress_01 = 0.;
		double sum_vm_plastic_strain = 0.;
		double sum_vm_stress = 0.;
		double sum_potential_energy = 0.;

		double sum2_stress_00 = 0.;
		double sum2_stress_11 = 0.;
		double sum2_stress_01 = 0.;
		double sum2_vm_plastic_strain = 0.;
		double sum2_vm_stress = 0.;
		double sum2_potential_energy = 0.;

		for(auto &element : elements)
		{
			auto &stress = element->stress();

			sum_stress_00 += stress[0][0];
			sum2_stress_00 += stress[0][0]*stress[0][0];

			sum_stress_11 += stress[1][1];
			sum2_stress_11 += stress[1][1]*stress[1][1];

			sum_stress_01 += stress[0][1];
			sum2_stress_01 += stress[0][1]*stress[0][1];

			double vm_plastic_strain = element->integrated_vm_eigenstrain();
			sum_vm_plastic_strain += vm_plastic_strain;
			sum2_vm_plastic_strain += vm_plastic_strain * vm_plastic_strain;

			double vm_stress = mepls::utils::get_von_mises_equivalent_stress(element->stress());
			sum_vm_stress += vm_stress;
			sum2_vm_stress += vm_stress * vm_stress;

			double potential_energy = 0.5 * dealii::invert(element->C()) * element->stress() * element->stress();
			sum_potential_energy += potential_energy;
			sum2_potential_energy += potential_energy * potential_energy;
		}

		double N = double(elements.size());

		data.av_stress_00 = sum_stress_00/N;
		data.av_stress_11 = sum_stress_11/N;
		data.av_stress_01 = sum_stress_01/N;
		data.av_vm_stress = sum_vm_stress/N;
		data.av_vm_plastic_strain = sum_vm_plastic_strain/N;
		data.av_potential_energy = sum_potential_energy/N;

		data.std_stress_00 = std::sqrt(
			sum2_stress_00/N - data.av_stress_00 * data.av_stress_00);

		data.std_stress_11 = std::sqrt(
			sum2_stress_11/N - data.av_stress_11 * data.av_stress_11);

		data.std_stress_01 = std::sqrt(
			sum2_stress_01/N - data.av_stress_01 * data.av_stress_01);

		data.std_vm_plastic_strain = std::sqrt(
			sum2_vm_plastic_strain/N - data.av_vm_plastic_strain * data.av_vm_plastic_strain);

		data.std_vm_stress = std::sqrt(sum2_vm_stress/N - data.av_vm_stress * data.av_vm_stress);

		data.std_potential_energy = std::sqrt(
			sum2_potential_energy/N - data.av_potential_energy * data.av_potential_energy);

		data.time = macrostate["time"];
		data.total_strain = macrostate["total_strain"];
		data.ext_stress = macrostate["ext_stress"];

		data.index = index;

      macro_evolution.push_back(data);
   }


	unsigned int current_index()
	{
		return index;
	}

	std::vector<DrivingRow> driving;
	/*!< Vector to store the information contained in objects of the class
	 * \ref event::Driving reformated as \ref History.DrivingRow. */

	std::vector<PlasticRow> plastic;
	/*!< Vector to store the information contained in objects of the class
	 * \ref event::Plastic reformated as \ref History.PlasticRow. */

	std::vector<RenewSlipRow> renew;
	/*!< Vector to store the information contained in objects of the class
	 * \ref event::RenewSlip reformated as \ref History.RenewSlipRow. */

	std::vector<MacroSummaryRow> macro_evolution;

	std::string name;

	unsigned int index;
	/*!< This member keeps track of the index to be assigned to the next
	 * registered event. It acts as a global identifier of each
	 * event, independent of its type. Hence, it corresponds to the order
	 * in which each change in the system occurs. The index is updated
	 * one unit in \ref History.add() each time an event is added to
	 * \ref History. Sometimes, it can have the same value for different
	 * events if those events are assumed to occur simultaneously.
	 * For example, if different thresholds within the same element are
	 * simultaneously renewed after plastic deformation occurs in that element.
	 *
	 * \note In the output datasets, it corresponds to the column "index".
	 * It allows reconstructing the simulation history by merging different
	 * datasets in the right way. */

	bool closed;
	/*!< If true, no further events are added. */
};


} // namespace history


} // namespace mepls

#endif //__HISTORY_H_
