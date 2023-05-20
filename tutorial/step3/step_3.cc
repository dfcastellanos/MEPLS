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


#include <example.h>
#include <mepls/utils.h>
#include <mepls/solver.h>
#include <mepls/system.h>
#include <mepls/event.h>

// new headers
#include <mepls/history.h>
#include <fstream>

int main()
{
	constexpr unsigned dim = 2;
	std::mt19937 generator(1234567);

	// let's consider a material with a shear modulus of 30 GPa,
	double G = 30.;
	double nu = 0.3;
	dealii::SymmetricTensor<4, dim> C = mepls::utils::tensor::make_isotropic_stiffness<dim>(G, nu);

	unsigned int Nx = 32;
	unsigned int Ny = 32;

	mepls::element::Vector<dim> elements;
	for(double n = 0; n < Nx * Ny; ++n)
	{
		example::element::Scalar<dim>::Config conf;
		conf.number = n;

		// plastic increments induce a local shear deformation of 5%
		conf.gamma = 0.05;

        // the slip thresholds have a scale of 1 GPa (the average is approx. the same, see Weibull distribution)
		conf.lambda = 1;

		// the slip thresholds have some disorder. For that, we use a not too-high k (for k->inf
		// the disorder vanishes, while it increases for k->0)
		conf.k = 6;

		auto element = new example::element::Scalar<dim>(conf, generator);
		element->C(C);

		elements.push_back(element);
	}

	mepls::elasticity_solver::LeesEdwards<dim> solver(Nx, Ny);
	for(auto &element : elements)
		solver.set_elastic_properties(element->number(), element->C());
	solver.setup_and_assembly();

	mepls::element::calculate_ext_stress_coefficients(elements, solver);
	mepls::element::calculate_local_stress_coefficients_central(elements, solver);

	mepls::system::Standard<dim> system(elements, solver, generator);


	// we create a history object to record the evolution of the sytem
	mepls::history::History<dim> sim_history("Simulation_history");

	// the system needs to know the history object so that it can record the driving
	// and slip events that we add
	system.set_history(sim_history);


	//----- DYNAMICS OF THE SYSTEM ------

	// record the system's initial macroscale properties
	sim_history.add_macro(system);

	// main simulation loop till 5% applied shear strain (epsilon_xy, not gamma)
	while(system.macrostate["total_strain"] < 0.05)
	{

		// print some output
		std::cout << system.macrostate["total_strain"] << " " << system.macrostate["ext_stress"]
				  << std::endl;

		//----- PERTURB THE SYTEM ------

		// the applied shear strain is (instantaneously) increased by a small
		// ammout of 0.01%, and is kept fixed afterward
		mepls::event::Driving<dim> driving_event;
		driving_event.dload = 0.0001;
		system.add(driving_event);

		// In this state the stress is at a peak, from where it will drop due to plastic
		// deformation
		sim_history.add_macro(system);


		//----- RELAX THE SYTEM ------

		std::vector<mepls::event::Plastic<dim> > events_relax_step;
		bool continue_relaxation = true;

		// relaxation loop (cascade of slip events)
		while(continue_relaxation)
		{
			events_relax_step.clear();

			// whenever we find an unstable slip system, we the corresponding plastic
			// event to the system. Instead of one at a time, we add a vector of them
			// to consider them simoultaneous in time (also, this increases performce since
			// the FEM problem is solved only once for all of them)
			for(auto &element : system)
				for(auto &slip : *element)
					if(slip->barrier < 0.)
					{
						mepls::event::Plastic<dim> plastic_event(slip);
						events_relax_step.push_back(plastic_event);

						// go to the next element, since the other slip system cannot have
						// a negative barrier too
						break;
					}

			system.add(events_relax_step);

			// if some slip events were added, the stress field has changed, and new slip systems
			// might be unstable, so we repeat the relaxation loop to check it. However, if no slip
			// event was added, which means the all the slip systems are stable, and no changes can occur
			// anymore. In this case, don't repeat the relaxation loop
			continue_relaxation = events_relax_step.size() > 0;

		} // relaxation loop

		// This state corresponds to a stress valley, after plastic deformation has stopped
		sim_history.add_macro(system);

	} // main simulation loop

	// remember that we allocated dynamically using the new operator
	// when we build a model, it will be our responsibility to delete them
	for(auto &element : elements)
		delete element;
}