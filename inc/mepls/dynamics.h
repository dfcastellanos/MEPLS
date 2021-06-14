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

#ifndef _DYNAMICS_H
#define _DYNAMICS_H

#include <mepls/system.h>
#include <mepls/snapshot.h>
#include <mepls/solver.h>

namespace mepls
{

/*! This namespace contains the algorithms to control the time-evolution of the
 * external load and the activation of slip systems within the material. */
namespace dynamics
{

/*! This namespace contains the values that indicate how an event
 * (see \ref event) was created (see \ref event::Plastic.activation_protocol
 * and \ref event::Driving.activation_protocol). */
namespace Protocol
{

enum Values
{
	variation_plasticity = 0, /*!< A driving event signaling an external stress or a total strain variation
    * induced by plastic activity. */

	ad_hoc = 1, /*!< An event explicitly created by the user. */

	extremal_dynamics = 2, /*!< A event created by \ref dynamics::extremal_dynamics_step. */

	kinetic_monte_carlo = 3, /*!< An event created by \ref dynamics::KMC. */

	relaxation_step = 4, /*!< An event created by \ref dynamics::relaxation. */

	prestress = 5,
	/*!< An driving event explicitly created by the user to represent the prestress. */

	variation_stiffness = 6, /*!< A driving event signaling a global stiffness variation. */

	metropolis_hastings = 7  /*!< An event created by \ref dynamics::MetropolisHastings. */
};

} // protocol


template<int dim>
std::vector<event::Plastic<dim>> relaxation(
								system::System<dim> &system,
								utils::ContinueSimulation &continue_simulation,
								double fracture_limit = 10);

/*! This class implements the Kinetic Monte Carlo method, used for simulating
 * the thermal activation of slip systems. */
template<int dim>
class KMC
{
  public:

	KMC()
		:
		unif_distribution(0, 1)
	{
		/*! Constructor. */
	}

	void operator()(system::System<dim> &system)
	{
		/*! Select a thermally activated slip system using the Kinetic Monte Carlo
		 * method and perform its associated plastic event. The selection is made
		 * using the activation rates of each slip system,
		 * \ref slip::Slip<dim>::activation_rate_. */

		/* -------- make a vector of rates and slips systems -------- */

		// we don't want to re-allocate the vectors each time. To this end, we
		// count the total number of slips systems first, and then call resize()
		// on the vectors. If the number of slip systems did not change, the
		// vectors are not re-allocated
		unsigned int n_slips = 0;
		for(auto &element : system)
			for(auto &slip : *element)
				++n_slips;

		double rate_sum = 0.;
		rates.resize(n_slips);
		slips.resize(n_slips);
		unsigned int n = 0;
		for(auto &element : system)
			for(auto &slip : *element)
			{
				rate_sum += slip->activation_rate;
				rates[n] = slip->activation_rate;
				slips[n] = slip;
				++n;
			}

		/* -------- compute time increment -------- */
		event::Driving<dim> load_increment_event;
		load_increment_event.dtime = utils::rand::get_exponential_rand(1. / rate_sum,
																	   unif_distribution(
																		   system.generator));
		load_increment_event.activation_protocol = Protocol::kinetic_monte_carlo;
		system.add(load_increment_event);


		/* -------- select the active slip system -------- */
		slip::Slip<dim> *thermal_slip;
		double rate_sum_old = 0.;
		double rate_sum_new = 0.;
		double r = unif_distribution(
			system.generator) * rate_sum; // unif. random between 0 and rate_sum
		for(unsigned int i = 0; i < slips.size(); ++i)
		{
			rate_sum_new += rates[i] + rate_sum_old;
			if(rate_sum_new > r)
			{
				// in this way, a slip system is chosen with a probability
				// proportional to its rate
				thermal_slip = slips[i];
				break;
			}
		}

		event::Plastic<dim> thermal_slip_event(thermal_slip);
		thermal_slip_event.activation_protocol = Protocol::kinetic_monte_carlo;
		system.add(thermal_slip_event);
	}

  private:
	std::uniform_real_distribution<double> unif_distribution;

	std::vector<mepls::slip::Slip<dim> *> slips;
	/*!< Vector containing all the slip systems within the material. It
	 * conatenates pointers to all the slip systems
	 * \ref element::Element<dim>::slip_systems_ in all the elements in \ref
	 * system.elements. */

	std::vector<double> rates;
	/*!< Vector containing the activation rates all the slip systems within the
	 * material. */
};





/*! This class implements the Metropolis-Hastings algorithm. */
template<int dim>
class MetropolisHastings
{
  public:

	MetropolisHastings()
		:
		unif_distribution(0, 1)
	{
		/*! Constructor. */
	}

	bool operator()(system::System<dim> &system, bool stress_redist = true)
	{
		/*! Perform as Metropolis-Hastings step: a random slip event is
		 * performed, and the global energy change computed. The event will be
		 * accepted if it energy change is negative. If it is positive, it will
		 * be accepted with probability exp(-energy_change/T). If stress_redist
		 * is false, the yielding element structural properties will be renewed
		 * but the elastic fields will not be updated.
		 *
		 * @return whether the event has been accepted or rejected. If it has
		 * been rejected, the state of the system is exactly the same as
		 * before calling this function. */

		auto & generator = system.generator;

		// we don't want to re-allocate the vectors each time. To this end, we
		// count the total number of slips systems first, and then call resize()
		// on the vectors. If the number of slip systems did not change, the
		// vectors are not re-allocated
		unsigned int n_slips = 0;
		for(auto &element : system)
			for(auto &slip : *element)
				++n_slips;

		slips.resize(n_slips);
		unsigned int n = 0;
		for(auto &element : system)
			for(auto &slip : *element)
			{
				slips[n] = slip;
				++n;
			}


		double energy_before = get_sum_energy(system.elements);

		auto stress_before = system.solver.get_stress();

		mepls::element::Vector<dim> elements_before;
		for(auto & element : system)
			elements_before.push_back( element->make_copy() );

		std::vector<event::Plastic<dim>> all_events;
		utils::ContinueSimulation continue_simulation;


		// select the candidate slip
		unsigned int i = int(unif_distribution(generator) * n_slips);
		auto candidate_slip = slips[i];
		event::Plastic<dim> cadidate_event( candidate_slip );
		cadidate_event.activation_protocol = Protocol::metropolis_hastings;
		auto element = candidate_slip->parent;

		if(stress_redist)
		{
			system.add(cadidate_event);
			all_events = dynamics::relaxation(system, continue_simulation);
			all_events.push_back(cadidate_event); // add also the first event
			if(not continue_simulation())
			{
				std::cout << continue_simulation << std::endl;
				abort();
			}
		}
		else
			element->renew_structural_properties();

		double energy_after = get_sum_energy(system.elements);

		// decide if it is accepted
		double energy_change = energy_after - energy_before;
		double prob_accept = energy_change < 0. ? 1. : std::exp(-energy_change/T);
		bool accepted = unif_distribution(generator) < prob_accept;

		// TODO if rejected and we called system.add, the last renew event
		//  history should be deleted

		if(not accepted)
		{
			for(auto & element : system)
				element->make_copy( elements_before[element->number()] );

			if(stress_redist)
			{
				// add the opposite eigenstrain increment to the solver, so that
				// in the next elastic solution the effects of this event are
				// undone. To keep things as clean as possible, do it also in the
				// element
				for(auto &event : all_events)
				{
					unsigned int n = event.element;
					system.solver.add_eigenstrain(n, -1.*event.eigenstrain);
				}

				// we avoid callig system.solve_elastic_problem() by setting the
				// previous	stress
				for(unsigned int n = 0; n < system.size(); ++n)
					system[n]->elastic_stress( stress_before[n] );
			}
		}

		for(auto &element : elements_before)
			delete element;

		return accepted;
	}


	double get_sum_energy(element::Vector<dim> &elements)
	{
		/*! Compute the average elastic energy over the input elements. */

		double sum_energy = 0.;

		for(auto &element : elements)
			sum_energy += element->energy();

		return sum_energy;
	}

	double T = 1.;

  private:
	std::uniform_real_distribution<double> unif_distribution;

	std::vector<mepls::slip::Slip<dim> *> slips;
	/*!< Vector containing all the slip systems within the material. It
	 * conatenates pointers to all the slip systems
	 * \ref element::Element<dim>::slip_systems_ in all the elements in \ref
	 * system.elements. */
};


template<int dim>
std::pair<event::Driving<dim>, event::Plastic<dim>> find_weakest_slip(
	const element::Vector<dim> &elements, bool is_forward)
{
	/*! Find the slip system whose activation (i.e., making
	 * slip::Slip<dim>::barrier < 0) requires the minimum increment of the load.
	 * This slip system is the least mechanically stable one in the material.
	 * If is_forward=True, the load increment is positive, otherwise it is
	 * negative.
	 *
	 * @return Returns a pair of events containing the change in the driving
	 * conditions (with the calculated minimum load increment) and the plastic
	 * event resulting from the activation of the slip system. */

	// sign of the load variation that we apply in order to unstabilize the slip
	// system
	double sign = is_forward ? 1. : -1.;

	mepls::slip::Slip<dim> *weakest_slip;
	event::Driving<dim> load_increment_event;
	load_increment_event.dload = sign * 1e90;
	double load_increment = 0.;

	for(auto &element : elements)
		for(auto &slip : *element)
		{
			load_increment = slip->get_critical_load_increment();
			if(std::abs(load_increment) < std::abs(
				load_increment_event.dload) and sign * load_increment > 0)
			{
				// Find which slip system requieres the smallest external load
				// increment to become unstable. The condition
				// sign*load_increment > 0 ensures that we are considering only
				// slip systems that require a load increment with the same sign
				// we are applying
				load_increment_event.dload = load_increment;
				weakest_slip = slip;
			}
		}

	M_Assert(not std::isnan(load_increment_event.dload), "");M_Assert(
		not(std::abs(load_increment_event.dload) > 1e20 and elements.size() == 1),
		"Maybe no positive load increment was found for a single-element-patch test?");M_Assert(
		std::abs(load_increment_event.dload) < 1e20,
		"Maybe no load increment with the right sign was found?");M_Assert(
		load_increment_event.dload != 0., "");M_Assert(sign * load_increment_event.dload > 0.,
													   "");M_Assert(weakest_slip != nullptr,
																	"Maybe no positive load increment was found?");M_Assert(
		weakest_slip->parent->number() < elements.size(), "");

	event::Plastic<dim> weakest_slip_event(weakest_slip);

	return std::make_pair(load_increment_event, weakest_slip_event);
}


template<int dim>
void extremal_dynamics_step(system::System<dim> &system, bool is_forward = true)
{
	/*! Perform the minimum external load increment necessary to activate one and
	 * only one slip system among all the existing ones in \ref system.elements.
	 * If is_forward=True, the load increment is positive; otherwise, it is
	 * negative.
	 *
	 * @note This function only performs the critical load increment but does not
	 * perform any plastic event. */

	auto ext_dyn_pair = find_weakest_slip(system.elements, is_forward);

	double sign = is_forward ? 1. : -1.;

	event::Driving<dim> &load_increment_event = ext_dyn_pair.first;

	// the computed load increment will make the slip system's barrier 0. We add
	// a small extra increment to make sure that the slip system will be
	// considered as (marginally) unstable by other algorithms (i.e. barrier < 0
	// although abs (barrier) is approx. 0 )
	load_increment_event.dload += sign * 1e-10;

	load_increment_event.activation_protocol = Protocol::extremal_dynamics;
	system.add(load_increment_event);

	// we don't trigger the plastic event here, only make the system marginally
	// stable

#ifdef DEBUG
	double min_barrier = 1e10;
	unsigned int n_neg_barriers = 0;
	for(auto &element : system)
	   for(auto &slip : *element)
	   {
		  if( slip->barrier < min_barrier )
			 min_barrier = slip->barrier;

		  if( slip->barrier < 0. )
			 n_neg_barriers += 1;

	   }
	M_Assert(sign*load_increment_event.dload > 0, "The load increment has the wrong sign");
	M_Assert(std::abs(load_increment_event.dload)<1e90, "The load increment value is infinite");
	M_Assert(min_barrier<0., "The minimum barrier after ext. dynamics is not negative");
	M_Assert(std::abs(min_barrier)<1e-6, "The minimum barrier after ext. dynamics is not small enough");
	M_Assert(n_neg_barriers==1., "Number of negative barriers after ext. dynamics is not equal to one");
#endif
}


template<int dim>
void finite_extremal_dynamics_step(
	double dicrete_increment, system::System<dim> &system, bool is_forward = true)
{
	/*! Perform external load increments in discrete steps of magnitude equal to
	 * the input dicrete_increment. Load increments are repeated until the first
	 * slip system among all the existing ones in \ref system.elements become
	 * active. Since the increments have a predetermined fixed value, at the
	 * moment of the first slip activation, more than one slip system might
	 * become simultaneously active, If is_forward=True, the load increment is
	 * positive; otherwise, it is negative.
	 *
	 * @note This function only performs the critical load increment but does not
	 * perform any plastic event. */


	auto ext_dyn_pair = dynamics::find_weakest_slip(system.elements, is_forward);
	event::Driving<dim> &load_increment_event = ext_dyn_pair.first;

	double sign = is_forward ? 1. : -1.;
	dicrete_increment *= sign;

	load_increment_event.activation_protocol = Protocol::extremal_dynamics;
	unsigned int n = std::floor(load_increment_event.dload / dicrete_increment);
	load_increment_event.dload = dicrete_increment * (n + 1);

	system.add(load_increment_event);
}


template<int dim>
void fixed_load_increment(double increment, system::System<dim> &system)
{
	/*! Perform an external load increment with value equal to the input
	 * dicrete_increment.
	 *
	 * @note This function only performs the critical load increment but does not
	 * perform any plastic event. */

	event::Driving<dim> load_increment_event;
	load_increment_event.activation_protocol = Protocol::ad_hoc;
	load_increment_event.dload = increment;

	system.add(load_increment_event);
}


template<int dim>
std::vector<event::Plastic<dim>> relaxation(
								system::System<dim> &system,
								utils::ContinueSimulation &continue_simulation,
								double fracture_limit)
{
	/*! Perform plastic events associated with each active slip system. The
	 * events are performed and the state of the system are updated. After that,
	 * new slip systems might become active as a consequence of changes in the
	 * elastic fields. If there is at least an active slip, the process is
	 * repeated. If not, the system is mechanically stable, and the function
	 * returns. */

	bool continue_relaxation = true;
	double ongoing_size = 0.;
	std::vector<event::Plastic<dim>> all_events;
	std::vector<event::Plastic<dim>> events_relax_step;

	/* ---- while there are mechanicaally unstable slip systems ----- */
	while(continue_relaxation and continue_simulation())
	{
		events_relax_step.clear();
		for(auto &element : system)
		{
			/* ---- find the weakest slip per element ----- */
			double lowest_barrier = 0.;
			mepls::slip::Slip<dim> *weakest_slip;
			for(auto &slip : *element)
				if(slip->barrier < lowest_barrier)
				{
					weakest_slip = slip;
					lowest_barrier = weakest_slip->barrier;
				}

			/* ---- activate the slip system if it is unstable ----- */
			if(lowest_barrier < 0.)
			{ M_Assert(weakest_slip != nullptr, "");
				event::Plastic<dim> plastic_event(weakest_slip);
				plastic_event.activation_protocol = Protocol::relaxation_step;
				events_relax_step.push_back(plastic_event);
			}
		}

		continue_relaxation = events_relax_step.size() > 0;
		ongoing_size += events_relax_step.size();

		/* ---- add & perform the events ----- */
		system.add(events_relax_step);

		for(auto & event : events_relax_step)
			all_events.push_back(event);


		continue_simulation(ongoing_size / double(system.size()) < fracture_limit,
							" fractured (avalanche size limit reached)");
	}

#ifdef DEBUG
	if(continue_simulation())
	   for(auto &element : system)
		  for(auto &slip : *element)
			 M_Assert(slip->barrier > 0, "Relaxation did not make all the barriers positive");
#endif

	return all_events;
}


} // namespace dynamcis
} // namespace discrete_platicity

#endif //_DYNAMICS_H
