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

#ifndef _SYSTEM_H
#define _SYSTEM_H

#include <random>
#include <math.h>

#include <mepls/solver.h>
#include <mepls/element.h>
#include <mepls/history.h>
#include <mepls/utils.h>

#include <deal.II/base/symmetric_tensor.h>

namespace mepls
{

namespace element
{
template<int dim>
class Element;
}


/*! This class represents the macroscopic state of the system. Its state is
 * updated by registering the occurence of external load variation events
 * (event::Driving) and plastic activity in the form of slip events
 * (event::Plastic). The values of the macroscopic magnitudes (with keys
 * "av_vm_plastic_strain", "load", "ext_stress", "pressure", "time" and
 * "total_strain") are then retrieved using the \ref operator[] with the key of
 * the magnitude. */
template<int dim>
struct MacroState
{
	MacroState(unsigned int n_elements_)
		:
		av_vm_plastic_strain(0.),
		total_strain(0.),
		load(0.),
		ext_stress(0.),
		pressure(0.),
		time(0.),
		n_elements(n_elements_)
	{
		/*! Constructor. */

		// Initialize the map of macroscopic magnitudes. The key-names given here
		// will be used in other parts of the code to query the state of the
		// different macroscale magnitudes
		monitor_map["av_vm_plastic_strain"] = &av_vm_plastic_strain;
		monitor_map["load"] = &load;
		monitor_map["ext_stress"] = &ext_stress;
		monitor_map["pressure"] = &pressure;
		monitor_map["time"] = &time;
		monitor_map["total_strain"] = &total_strain;
	};

	MacroState & operator=(const MacroState &rhs)
	{
		/*! Assigment operator. Copy the rhs object into the lhs. */

		av_vm_plastic_strain = rhs.av_vm_plastic_strain;
		total_strain = rhs.total_strain;
		load = rhs.load;
		ext_stress = rhs.ext_stress;
		pressure = rhs.pressure;
		time = rhs.time;
		n_elements = rhs.n_elements;

		// monitor_map must not be copied, otherwise its pointers would
		// point to the members of the object from where the copy is made

		return *this;
	};

	MacroState(const MacroState &input_macrostate)
	{
		/*! Copy constructor. */

		// use the constructor to initialize the monitor_map
		MacroState( input_macrostate.n_elements );

		// use the = operator to copy the state of all the members except
		// the monitor map
		*this = input_macrostate;
	};

	double operator[](const std::string &name) const
	{
		/*! Return the current value of the macroscopic magnitude with the given
		 * name. */

		// The input key-name is used as key for \ref monitor_map. Since \ref
		// monitor_map stores pointers, the returned value is the current
		// up-to-date value.
		M_Assert(monitor_map.find(name) != monitor_map.end(), "Monitor magnitude does not exist");

		return *monitor_map.at(name);
	};

	void update(event::Driving<dim> &event)
	{
		/*! Update the macroscopic quantities using the information contained in
		 * the input \ref event::Driving object. */

		load += event.dload;
		ext_stress += event.dext_stress;
		pressure += event.dpressure;
		total_strain += event.dtotal_strain;
		time += event.dtime;
	};

	void update(event::Plastic<dim> &event)
	{
		/*! Update the macroscopic quantities using the information contained in
		 * the input \ref event::Plastic object. */

		av_vm_plastic_strain += event.dplastic_strain / double(n_elements);
		time += event.dtime;
	};

	void clear()
	{
		/*! Reset the value of all the magnitudes to zero. */

		av_vm_plastic_strain = 0.;
		total_strain = 0.;
		load = 0.;
		ext_stress = 0.;
		time = 0.;
	};

  private:
	double av_vm_plastic_strain;
	/*!< Average plastic strain. It is defined as the sum of the von Mises
	 * equivalent plastic strain introduced by all events slip events, divided by
	 * the number of mesoscale elements in the system.
	 *
	 * @note it is not the net plastic deformation, i.e., adding events of
	 * magnitude 1 and -1 does not result in 0 average plastic strain but in
	 * +2. */

	double total_strain;
	/*!< Total strain (i.e. the sum of elastic and plasic strain) undergone by
	 * the system. Its relation with \ref load depends on the specific loading
	 * conditions (see \ref
	 * elasticity_solver::ElasticitySolver.get_total_strain()). */

	double load;
	/*!< Value of the applied load. It refers to the amount of imposed
	 * displacement (if displacement-controlled) or the modulus of the applied
	 * traction (if traction-controlled). The control mode is defined by \ref
	 * base_param::Simulation.control_mode. */

	double pressure;
	/*!< Average pressure in the system. */

	double ext_stress;
	/*!< Average of the relevant component of the external stress tensor along
	 * the loaded surface or surfaces. Its relation with \ref load and its
	 * computations depends on the loading mode (see \ref
	 * elasticity_solver::ElasticitySolver.get_loading_stress()). */

	double time;
	/*!< Elapsed time. */

	std::map<std::string, const double *> monitor_map;
	/*!< Map relating the names of the macroscopic magnitudes to pointers to
	 * their values. Consequently, when the value of a magnitude with a certain
	 * name is requested, the map returns the current value of the magnitude.
	 * Existing keys are "av_vm_plastic_strain", "load", "ext_stress", "pressure",
	 * "time" and "total_strain". */

	unsigned int n_elements;
	/*!< Number of elements composing the system. */

};



/*! This namespace contains the abstract class system::System, which
 * serves as the interface for the system objects. A system object represents a
 * solid material capable of undergoing plastic deformation, which occurs in
 * a discrete manner both in time and in space. */
namespace system
{


/*! The System class represents a solid material capable of undergoing plastic
 * deformation, which occurs in a discrete manner both in time and in space.
 * The class is composed of different members that represent different aspects
 * of the material: a vector \ref elements of element objects which represent
 * discrete non-overlapping mesoscale material regions; an elasticity solver
 * (\ref solver) to compute the elastic fields; the member \ref history to
 * register the details of each event occurring in the system. Normally,
 * system objects are operated by the algorithms in \ref dynamics controlling
 * its evolution using the events defined in \ref event. The operation is done
 * by calling the \ref add() functions. */
template<int dim>
class System
{
  public:

	System(
		element::Vector<dim> &elements_,
		elasticity_solver::Solver<dim> &solver_,
		std::mt19937 &generator_)
		:
		generator(generator_),
		macrostate(elements_.size()),
		history(&default_history),
		elements(elements_),
		solver(solver_)
	{
		/*! Constructor. Initialized reference members. */
	}

	virtual ~System()
	{
		/*! Virtual destructor. */
	}

	typename element::Vector<dim>::iterator begin()
	{
		/*! Beginning of the iteration over the system. It corresponds to the
		 * beginning of the vector containing the elements composing the system.*/

		return elements.begin();
	}

	typename element::Vector<dim>::iterator end()
	{
		/*! End of the iteration over the system. It corresponds to the end of the
		 * vector containing the elements composing in the system. */

		return elements.end();
	}

	typename element::Vector<dim>::iterator begin() const
	{
		/*! Beginning of the iteration over the system. It corresponds to the
		 * beginning of the vector containing the elements composing the system.*/

		return elements.begin();
	}

	typename element::Vector<dim>::iterator end() const
	{
		/*! End of the iteration over the system. It corresponds to the end of the
		 * vector containing the elements composing in the system. */

		return elements.end();
	}

	unsigned int size() const
	{
		return elements.size();
	}

	virtual void solve_elastic_problem(
		const std::vector<event::Plastic<dim>> &added_yielding,
		event::Driving<dim> &driving_event) = 0;
	/*!< This function defines how the elastic fields are computed and delivered
	 * to each element. */

	virtual void add(event::Driving<dim> &driving_event) = 0;
	/*!< This function defines how driving events are performed. */

	virtual void add(std::vector<event::Plastic<dim>> &added_yielding) = 0;

	/*!< This function defines how plastic (slip) events are performed. */

	void add(event::Plastic<dim> &plastic_event)
	{
		/*! This function wraps \ref add(std::vector<event::Plastic<dim>> &) to
		 * allow the user to add a single plastic event without the need of
		 * creating a vector.
		 *
		 * @note a vector with the input event is creted during the function call.
		 * If many events are to be added at once, \ref
		 * add(std::vector<event::Plastic<dim>> &) is the preferred way. */

		// TODO we should pass the original input plastic event instead of a copy,
		// so after this call the user can inspect the changes in the event
		// (i.e. the information about the event which is added during the call
		// to add)
		std::vector<event::Plastic<dim>> added_yielding = {plastic_event};
		add(added_yielding);
	}

	virtual System<dim> *get_new_instance(
		element::Vector<dim> &elements_,
		elasticity_solver::Solver<dim> &solver_,
		std::mt19937 &generator_) const = 0;

	/*!< This function creates a new object of the same dynamic type as the one
	 * from where it is called. The input parameters are used for the
	 * construction of the new object.
	 *
	 * @warning the user is responsible for deleting the returned object. */


	void set_history(history::History<dim> &history_)
	{
		history = &history_;
	}

	std::mt19937 &generator;
	/*!< Random number generator. */

	MacroState<dim> macrostate;
	/*!< Macrosopic state of the system. */

	history::History<dim> default_history;

	history::History<dim> *history;
	/*!< History of added plastic and driving events. */

	element::Vector<dim> &elements;
	/*!< Elements composing the system. */

	elasticity_solver::Solver<dim> &solver;
	/*!< Elasticity solver used for computing the elastic fields. */
};


/*! This class implements the most basic and widely-applicable functionality for
 * the system. See \ref solve_elastic_problem() and \ref add() functions for
 * details. */
template<int dim>
class Standard: public System<dim>
{

  public:

	Standard(
		element::Vector<dim> &elements_,
		elasticity_solver::Solver<dim> &solver_,
		std::mt19937 &generator_)
		:
		System<dim>(elements_, solver_, generator_)
	{
		/*! Constructor. */
	}

	system::System<dim> *get_new_instance(
		element::Vector<dim> &elements_,
		elasticity_solver::Solver<dim> &solver_,
		std::mt19937 &generator_) const override
	{
		return new Standard<dim>(elements_, solver_, generator_);
	}


	void solve_elastic_problem(
		const std::vector<event::Plastic<dim>> &added_yielding,
		event::Driving<dim> &driving_event) override
	{
		/*! This function uses the \ref solver to compute the elastic equilibrium
		 * state using for the eigenstrain field stored in the elements under the
		 * action of the current external load value. The stress field of each
		 * element is updated accordingly. The macroscale properties stored in
		 * \ref macrostate are also updated.
		 *
		 * While the value of the external load is established outside the system,
		 * other external conditions can suffer changes by getting feedback from
		 * the operations performed within the system. As an example of such
		 * variation, if the system is displacement-controlled, plastic
		 * deformation will cause the external stress to drop. On the other hand,
		 * if it is traction-controlled, it will cause the total strain to
		 * increase. Those changes are stored in the input driving event. */

		solver.solve();

		auto &stress = solver.get_stress();
		for(auto &element : elements)
		{
			unsigned int n = element->number();
			element->elastic_stress(stress[n]);

			// We set the stress back to the solver because the element adds to the
			// elastic stress the prestress. In this way, the solver will take into
			// account the total stress for computing the external stress instead
			// of jus the elastic stress. An alternative would be to let the solver
			// own the prestress tensor. However, the prestress is set, changed, and
			// reset to zero when manipulating the elements in different parts of
			// the code, so it is much easier to do it in this way.
			solver.set_stress(n, element->stress());
		}

		double external_stress = solver.get_external_stress();

		// compute the average pressure
		double pressure = 0.;
		for(auto &element : elements)
			pressure += -dealii::trace(element->stress()) / double(dim);
		pressure /= double(elements.size());

		// compute the changes in the driving mechanism with respect to the
		// previous elastic equilibrium state
		driving_event.dpressure = pressure - macrostate["pressure"];
		driving_event.dext_stress = external_stress - macrostate["ext_stress"];
		driving_event.dtotal_strain = solver.get_total_strain() - macrostate["total_strain"];
	}


	void add(event::Driving<dim> &driving_event) override
	{
		/*! Update the external load by adding the load increment defined in the
		 * input driving event (see event::Driving<dim>::dload). The elastic
		 * fields are computed with \ref solver, and the \ref elements are informed
		 * about their values. The \ref macrostate is updated. The input driving
		 * event is registered in the \ref history. */

		solver.add_load_increment(driving_event.dload);

		std::vector<event::Plastic<dim>> added_yielding;

		solve_elastic_problem(added_yielding, driving_event);

		history->add(driving_event);

		macrostate.update(driving_event);
	}


	void add(std::vector<event::Plastic<dim>> &added_yielding) override
	{
		/*! Update the plastic (eigenstrain) field of the \ref elements using the
		 * input vector of plastic events. Each plastic event defines the
		 * activation of a specific slip system (see \ref slip). For all the
		 * active slip systems, an eigenstrain increment is computed with \ref
		 * slip::Slip<dim>::get_eigenstrain_increment(), which is then added to
		 * their respective parent elements (see \ref slip::Slip<dim>::parent).
		 * The elastic fields are computed with \ref solver and the \ref elements
		 * are informed about their values. The \ref macrostate is updated. The
		 * input plastic events are registered in the \ref history. All the
		 * elements in which a slip event took place get their structural
		 * properties renewed, as defined by \ref
		 * element::Element<dim>::renew_structural_properties(). The new structural
		 * properties are resgitered as events of type \ref
		 * event::RenewSlip<dim> in \ref history.
		 *
		 * @note since plastic deformation changes the external stress (if the
		 * system is displacement-controlled) or the total strain (if the system
		 * is traction-controlled), such changes in the external conditions are
		 * registered in the \ref history as driving events
		 * (\ref event::Driving).*/

		if(added_yielding.size() == 0)
			return;

		// Iterate over the plastic events and: calculate their equivalent
		// eigenstrain, update the element eigenstrain tensor, upate the
		// macroscopic state, update the elasticity solver, and put the plastic
		// events from the map into a vector so they can be added simultaneously
		// to the history with the same output index
		std::vector<event::Plastic<dim>> combined_platic_events;
		std::vector<element::RenewInstruct<dim>> combined_renew_instruct;
		for(auto &plastic_event : added_yielding)
		{
			auto *slip = plastic_event.slip;
			auto *element = slip->parent;

			dealii::SymmetricTensor<2, dim> eigenstrain_increment = slip->get_eigenstrain_increment();

			plastic_event.eigenstrain = eigenstrain_increment;
			plastic_event.dplastic_strain = utils::get_von_mises_equivalent_strain(
				eigenstrain_increment);

			element::RenewInstruct<dim> renew_instruct;
			renew_instruct.slip_properties = true;
			renew_instruct.elastic_properties = false;
			renew_instruct.plastic_event = plastic_event;

			element->add_eigenstrain(eigenstrain_increment);

			macrostate.update(plastic_event);

			solver.add_eigenstrain(element->number(), eigenstrain_increment);

			combined_renew_instruct.push_back(renew_instruct);
			combined_platic_events.push_back(plastic_event);
		}
		history->add(combined_platic_events);



		// update the variation of the external driving due to platic deformation
		// (external stress drops occur if the system is displacement-controlled)
		event::Driving<dim> driving_variation_event;
		solve_elastic_problem(added_yielding, driving_variation_event);
		history->add(driving_variation_event);
		macrostate.update(driving_variation_event);

#ifdef DEBUG
		// sometimes external stress increments occur after plastic deformation.
		// This can sometimes happen if the system is very small. However, big
		// increments mean that there is a bug.
		if(driving_variation_event.dext_stress > 1e-10)
		   std::cout << "* Warning: plastic events made the external stress rise by "
					 << driving_variation_event.dext_stress << std::endl;
#endif


		// renew the structural properties of the elements in which a plastic
		// event has taken place
		std::vector<event::RenewSlip<dim>> renewal_vector;
		for(auto &r : combined_renew_instruct)
		{
			auto *slip = r.plastic_event.slip;
			auto *element = slip->parent;
			element->renew_structural_properties(r);
			element->record_structural_properties(renewal_vector);
		}
		history->add(renewal_vector);

		added_yielding.clear();
	}

	using system::System<dim>::macrostate;
	using system::System<dim>::history;
	using system::System<dim>::generator;
	using system::System<dim>::elements;
	using system::System<dim>::solver;
	using system::System<dim>::solve_elastic_problem;
	using system::System<dim>::add;

};


/*! This class implements the same functionality as \ref Standard<dim> but
 * considers randomized elastic interactions between the \ref elements. This is
 * achieved by random-shuffling the stress variations of the \ref elements when
 * their elastic fields are updated. In this way, we remove spatial
 * correlations in the variations but retain other statistical properties. */
template<int dim>
class ShuffledKernel: public Standard<dim>
{

  public:

	ShuffledKernel(
		element::Vector<dim> &elements_,
		elasticity_solver::Solver<dim> &solver_,
		std::mt19937 &generator_)
		:
		Standard<dim>(elements_, solver_, generator_)
	{
		/*! Constructor. */
	}

	system::System<dim> *get_new_instance(
		element::Vector<dim> &elements_,
		elasticity_solver::Solver<dim> &solver_,
		std::mt19937 &generator_) const override
	{
		return new ShuffledKernel<dim>(elements_, solver_, generator_);
	}


	void solve_elastic_problem(
		const std::vector<event::Plastic<dim>> &added_yielding,
		event::Driving<dim> &driving_event) override
	{
		/*! Perform the same operations as Standard<dim>::solve_elastic_problem.
		 * However, the \ref elements interact via randomized elastic fields. To
		 * this end, the changes in the elastic fields are computed and then
		 * element-wise random-shuffled. In this way, the changes retain all their
		 * statistical properties except for spatial correlations.
		 *
		 * @note The changes in the elastic fields are the consequence of slip
		 * events that occur within certain \ref elements. Elements undergoing slip
		 * events are not considered in the random-shuffling of elastic fields to
		 * ensure that slip events do reduce local stress fields in the mesoscale
		 * elements in which they occur.
		 *
		 * @warning The random-shuffling affects the total local change, that is,
		 * the superposition of external and internal stress fields. Thus,
		 * the shuffling operation is equivalent to a randomized interaction
		 * kernel only if the loading mode induces a homogeneous external stress
		 * field (as, e.g., in the case of pure shear or tension). */

		solver.solve();

		// since we call \ref solver.clear() at the end of this function, the call
		// to \ref solver.get_stress() returns the stress change instead of the
		// total stress
		auto &stress_change = solver.get_stress();

		// find the numbers of the yielded elements (i.e. that undergo plastic
		// deformation)
		std::set<unsigned int> yielded_elements;
		for(auto &plastic_event : added_yielding)
			yielded_elements.insert(plastic_event.slip->parent->number());

		// find the numbers of the elements that did not undergo plastic
		// deformation
		std::vector<unsigned int> non_yielded_elements;
		for(auto &element : elements)
		{
			unsigned int n = element->number();
			if(not yielded_elements.count(n))
				non_yielded_elements.push_back(n);
		}

		// shuffle the order of the numbers of the non-yielded elements.
		std::shuffle(non_yielded_elements.begin(), non_yielded_elements.end(), generator);


		unsigned int nn = 0;
		for(auto &element : elements)
		{
			unsigned int n = element->number();

			// If the element number #n did't yield, the stress is updated by
			// adding to its current value the change of the element number #m,
			// where m = non_yielded_elements[nn]. If the vector
			// non_yielded_elements was not shuffled then n = m, but in this case
			// n != m, fulfilling that element m did not yield). If the element
			// number #n did yield, its stress is updated normally.
			if(not yielded_elements.count(n))
				element->elastic_stress(
					element->elastic_stress() + stress_change[non_yielded_elements[nn++]]);
			else
				element->elastic_stress(element->elastic_stress() + stress_change[n]);

			// We set the stress back to the solver because the element adds to the
			// elastic stress the prestress. In this way, the solver will take into
			// account the total stress for computing the external stress instead
			// of just the elastic stress. An alternative would be to let the
			// solver own the prestress tensor. However, the prestress is set,
			// changed, and reset to zero when manipulating the elements in
			// different parts of the code, so it is much easier to do it in this
			// way.
			solver.set_stress(n, element->stress());
		}


		// update macroscate properties in the same way as in \ref
		// system::System<dim>::solve_elastic_problem()
		double external_stress = solver.get_external_stress();

		double pressure = 0.;
		auto stress_tensors = solver.get_stress();
		for(auto &tensor : stress_tensors)
			pressure += dealii::trace(tensor) / double(dim);
		pressure /= double(stress_tensors.size());

		driving_event.dpressure = pressure - macrostate["pressure"];
		driving_event.dext_stress = external_stress - macrostate["ext_stress"];

		// since we call \ref solver.clear() at the end of this function, we
		// cannot ask the \ref solver for the total strain. Under
		// displacement-controlled conditions, the external strain increment is
		// defined in the input driving_event.dload.

		// TODO the traction-controlled case should be considered. In that case,
		// is driving_event.dext_stress what is given by driving_event.dload.
		// A way to achieve this would be to let system::System<dim> know the
		// loading control mode, or to as \ref solver about it.
		driving_event.dtotal_strain = driving_event.dload;


		// clear the state of \ref solver.clear() so that when we call \ref
		// solver.get_stress() at the beginning of this function we obtain the
		// stress change isntead of the total stress
		solver.clear();
	}

	using system::Standard<dim>::macrostate;
	using system::Standard<dim>::history;
	using system::Standard<dim>::generator;
	using system::Standard<dim>::elements;
	using system::Standard<dim>::solver;
};


/*! This class implements the same functionality as \ref Standard<dim> but
 * considers homgeneous elastic interactions between the \ref elements
 * (mean-field interaction). This is achived by computing the average stress
 * variations over the \ref elements when their elastic fields are updated. */
template<int dim>
class HomogeneousKernel: public Standard<dim>
{

  public:

	HomogeneousKernel(
		element::Vector<dim> &elements_,
		elasticity_solver::Solver<dim> &solver_,
		std::mt19937 &generator_)
		:
		Standard<dim>(elements_, solver_, generator_)
	{
		/*! Constructor. */
	}

	system::System<dim> *get_new_instance(
		element::Vector<dim> &elements_,
		elasticity_solver::Solver<dim> &solver_,
		std::mt19937 &generator_) const override
	{
		return new HomogeneousKernel<dim>(elements_, solver_, generator_);
	}


	void solve_elastic_problem(
		const std::vector<event::Plastic<dim>> &added_yielding,
		event::Driving<dim> &driving_event) override
	{
		/*! Perform the same operations as Standard<dim>::solve_elastic_problem.
		 * However, the \ref elements interact via homogeneous elastic fields. To
		 * this end, the changes in the elastic fields are computed and then
		 * averaged.
		 *
		 * @note The changes in the elastic fields are the consequence of slip
		 * events that occur within certain \ref elements. Elements undergoing slip
		 * events are not considered in the random-shuffling of elastic fields to
		 * ensure that slip events do reduce local stress fields in the mesoscale
		 * elements in which they occur.
		 *
		 * @warning The averaging affects the total local change, that is, the
		 * superposition of external and internal stress fields. Thus, the
		 * averaging operation is equivalent to a homogeneous (mean-field)
		 * interaction kernel only if the loading mode induces a homogeneous
		 * external stress field (as, e.g., in the case of pure shear or
		 * tension). */

		solver.solve();

		// since we call \ref solver.clear() at the end of this function, the call
		// to \ref solver.get_stress() returns the stress change instead of the
		// total stress
		auto &elastic_stress_increment = solver.get_stress();

		// find the numbers of the yielded elements (i.e. that undergo plastic
		// deformation)
		std::set<unsigned int> yielded_elements;
		for(auto &plastic_event : added_yielding)
			yielded_elements.insert(plastic_event.slip->parent->number());

		// find the numbers of the elements that did not undergo plastic
		// deformation and compute the average stress change
		dealii::SymmetricTensor<2, dim> av_stress_increment;
		std::set<unsigned int> non_yielded_elements;
		for(auto &element : elements)
		{
			unsigned int n = element->number();
			if(not yielded_elements.count(n))
			{
				non_yielded_elements.insert(n);
				av_stress_increment += elastic_stress_increment[n];
			}
		}
		av_stress_increment /= double(non_yielded_elements.size());


		for(auto &element : elements)
		{
			unsigned int n = element->number();

			// If the element number #n did't yield, the stress is updated by
			// adding to its current value the change of the element number #m,
			// where m = non_yielded_elements[nn]. If the vector
			// non_yielded_elements was not shuffled then n = m, but in this case
			// n != m, fulfilling that element m did not yield). If the element
			// number #n did yield, its stress is updated normally.
			if(not yielded_elements.count(n))
				element->elastic_stress(element->elastic_stress() + av_stress_increment);
			else
				element->elastic_stress(element->elastic_stress() + elastic_stress_increment[n]);

			// We set the stress back to the solver because the element adds to
			// the elastic stress the prestress. In this way, the solver will take
			// into account the total stress for computing the external stress
			// instead of just the elastic stress. An alternative would be to let
			// the solver owns the prestress tensor. However, the prestress is set,
			// changed, and reset to zero when manipulating the elements in
			// different parts of the code, so it is much easier to do it in this
			// way.
			solver.set_stress(n, element->stress());
		}

		// update macroscate properties in the same way as in \ref
		// system::System<dim>::solve_elastic_problem()
		double external_stress = solver.get_external_stress();

		double pressure = 0.;
		auto stress_tensors = solver.get_stress();
		for(auto &tensor : stress_tensors)
			pressure += dealii::trace(tensor) / double(dim);
		pressure /= double(stress_tensors.size());

		driving_event.dpressure = pressure - macrostate["pressure"];
		driving_event.dext_stress = external_stress - macrostate["ext_stress"];

		// since we call \ref solver.clear() at the end of this function, we
		// cannot ask the \ref solver for the total strain. Under
		// displacement-controlled conditions, the external strain increment is
		// defined in the input driving_event.dload.

		// TODO the traction-controlled case should be considered. In that case,
		// is driving_event.dext_stress what is given by driving_event.dload. A
		// way to achieve this would be to let system::System<dim> know the
		// loading control mode, or to as \ref solver about it.
		driving_event.dtotal_strain = driving_event.dload;

		// clear the state of \ref solver.clear() so that when we call \ref
		// solver.get_stress() at the beginning of this function we obtain the
		// stress change instead of the total stress
		solver.clear();
	}

	using system::Standard<dim>::macrostate;
	using system::Standard<dim>::history;
	using system::Standard<dim>::generator;
	using system::Standard<dim>::elements;
	using system::Standard<dim>::solver;
};

} // namespace system

} // namespace mepls

#endif //_SYSTEM_H
