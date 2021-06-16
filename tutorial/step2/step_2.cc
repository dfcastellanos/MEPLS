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

// new headers
#include <mepls/solver.h>
#include <mepls/system.h>
#include <mepls/history.h>

int main()
{

	// ------ SETTING UP THE ELEMENTS ------

	constexpr unsigned dim = 2;
    std::mt19937 generator(1234567);

	double G = 1.;
	double nu = 0.3;
	dealii::SymmetricTensor<4,dim> C = mepls::utils::tensor::make_isotropic_stiffness<dim>(G, nu);

	unsigned int Nx = 16;
	unsigned int Ny = 16;

    // this type is an alias for std::vector<mepls::element::Element<dim>>
	mepls::element::Vector<dim> elements;

	for(double n = 0; n < Nx * Ny; ++n)
	{
		example::element::Scalar<dim>::Config conf;
		conf.number = n;

		auto element = new example::element::Scalar<dim>(conf, generator);
		element->C( C );

		elements.push_back(element);
	}

	// create the solver
	mepls::elasticity_solver::LeesEdwards<dim> solver(Nx, Ny);

	// set the elastic properties using the ones specified for the elements
	for(auto &element : elements)
		solver.set_elastic_properties(element->number(), element->C());

    // this operation builds the FEM system. After this is called, the solver is
    // ready to be used and its mesh geometry or loading model cannot be altered
	solver.setup_and_assembly();



	// ------ EXPLING THE SOLVER'S INTERFACE ------

	std::cout << "Testing the solver..." << std::endl;

	double delta_external_epsilon_xy = 0.1;
    solver.add_load_increment( delta_external_epsilon_xy );

    dealii::SymmetricTensor<2,dim> local_eigenstrain_increment;
	local_eigenstrain_increment[0][1] = 0.2;
	unsigned int element_number = 120;
	solver.add_eigenstrain(element_number, local_eigenstrain_increment);

	solver.solve();

	const std::vector<dealii::SymmetricTensor<2,dim>> & stress_field = solver.get_stress();
	const std::vector<dealii::Tensor<2,dim>> & def_grad = solver.get_deformation_gradient();

    double external_stress = solver.get_external_stress();
    double total_strain = solver.get_total_strain();

	// just tau = G * 2eps_xy
	std::cout << "    * external stress = " << external_stress << std::endl;

	// just eps_xy
	std::cout << "    * total strain = " << total_strain << std::endl;

	std::cout << "Cleaning the solver..." << std::endl;

	// this calls resets the solver state so the changes we made to load and the
	// eigenstrain field dissapear
	solver.clear();


	// ------ FINISHING THE ELEMENTS SETUP ------

	mepls::element::calculate_ext_stress_coefficients(elements, solver);

	// let's take some element
	auto element = elements[40];

	// let's take the first slip of that element
	auto slip = element->slip(0);

	// this is the external load increment that is necessary to make this slip's sytem
	// shear stress match its threshold, that is to make the barrier equal 0
	double critical_incr = slip->get_critical_load_increment();

	std::cout << "\nSlip's system critical load increment = " << critical_incr << std::endl;


	mepls::element::calculate_local_stress_coefficients_central(elements, solver);

    // this increment is to be added to the slips parent element to represents the effects
    // of the slip event
	dealii::SymmetricTensor<2,dim> eigenstrain_incr = slip->get_eigenstrain_increment();

	std::cout << "Slip's system local eigenstrain increment = " << eigenstrain_incr << std::endl;


	// ------ CREATING A SYSTEM ------

	mepls::system::Standard<dim> system(elements, solver, generator);

	auto & macrostate = system.macrostate;


	// ------ ADDING A SLIP EVENT ------

	mepls::event::Plastic<dim> plastic_event( slip );
	system.add(plastic_event);

	std::cout << "\nAfter adding the slip event: "
			  << "\n    * Local eigenstrain = " << element->eigenstrain()
			  << "\n    * Local von Mises eigenstrain = " << element->integrated_vm_eigenstrain()
			  << "\n    * Local stress = " << element->stress()
			  << "\n    * External stress = " << macrostate["ext_stress"]
			  << "\n    * Total strain = " << macrostate["total_strain"]
			  << "\n    * Global von Mises plastic strain = " << macrostate["av_vm_plastic_strain"]
			  << std::endl;

	system.solver.write_vtu("slip_event.vtu");


	// ------ ADDING A DRIVING EVENT ------

    mepls::event::Driving<dim> driving_event;

    // a load increment of 0.01
    driving_event.dload = 0.01;

    system.add(driving_event);


	std::cout << "\nAfter adding the driving event: "
			  << "\n    * Local eigenstrain = " << element->eigenstrain()
			  << "\n    * Local von Mises eigenstrain = " << element->integrated_vm_eigenstrain()
			  << "\n    * Local stress = " << element->stress()
			  << "\n    * External stress = " << macrostate["ext_stress"]
			  << "\n    * Total strain = " << macrostate["total_strain"]
			  << "\n    * Global von Mises plastic strain = " << macrostate["av_vm_plastic_strain"]
			  << std::endl;

	system.solver.write_vtu("slip_and_driving_events.vtu");


	// remember that we allocated dynamically using the new operator
	// when we build a model, it will be our responsibility to delete them
	for(auto &element : elements)
		delete element;
}