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

#include <deal.II/base/logstream.h>
#include <deal.II/base/path_search.h>

#include "example.h"

#include <cmdparser.hpp>


namespace example
{

using namespace mepls;


template<int dim>
class Simulation : public mepls::utils::Launcher<parameters::Parameters>
{

void run_impl(const parameters::Parameters &p) override
{

	/////////////////////////////
	//////  setup system ///////
	///////////////////////////


	auto timer = utils::TimerSingleton::getInstance();

	timer->enter_subsection("Setting up");

	std::mt19937 generator(p.sim.seed);

	auto C = utils::tensor::make_isotropic_stiffness<dim>(p.mat.shear_modulus,
	                                                       p.mat.poissons_ratio);

	std::vector<mepls::element::Element<dim> *> elements;
	for(double n = 0; n < p.sim.Nx * p.sim.Ny; ++n)
	{
		typename element::Example<dim>::Config conf;
		conf.threshold_shape = p.mat.threshold_shape;
		conf.threshold_scale_quench = p.mat.threshold_scale_quench;
		conf.threshold_scale = p.mat.threshold_scale;
		conf.coupling_constant = p.mat.coupling_constant;
		conf.activation_rate = p.mat.activation_rate;
		conf.C = C;
		conf.number = n;

		elements.push_back(
			new element::Example<dim>(generator, conf));
	}

	elasticity_solver::LeesEdwards<dim> solver(p.sim.Nx, p.sim.Ny);
	M_Assert(elements.size() == solver.get_n_elements(), "Numbers of mesoscale and finite "
													  "elements do not match");
	for(auto &element : elements)
		solver.set_elastic_properties(element->number(), element->C());
	solver.setup_and_assembly();

	mepls::element::calculate_local_stress_coefficients_central(elements, solver);
	mepls::element::calculate_ext_stress_coefficients(elements, solver);



	utils::ContinueSimulation continue_simulation;

	std::vector<snapshot::Stress<dim> > stress_snapshots;
	snapshot::Check<dim> snapshot_check(p.out.snapshots_min, p.out.snapshots_max,
										p.out.snapshots_interval, p.out.snapshots_sensitivity);
										
	system::Standard<dim> system(elements, solver, generator);
	MacroState<dim> &macrostate = system.macrostate;

	timer->leave_subsection("Setting up");



	/////////////////////////////////
	//////      dynamics     ///////
	///////////////////////////////

	history::History<dim> history("athermal_quasistatic_shear");
	system.set_history(history);

	/* ----- snapshots ----- */
	stress_snapshots.push_back(
		snapshot::Stress(system, "total_strain", 0.,
						 macrostate["total_strain"]));



	/* ---- simulation loop ----- */
	timer->enter_subsection("Dynamics");


	history.add_macro(system);

	while(continue_simulation())
	{

		if(p.out.verbosity and omp_get_thread_num() == 0)
			std::cout << history.index() << " | " << std::fixed << macrostate["total_strain"]
					  << " " << macrostate["ext_stress"] << " " << macrostate["pressure"]
					  << std::endl;

		if(snapshot_check(macrostate["total_strain"]))
		{

			stress_snapshots.push_back(
				snapshot::Stress(system, "total_strain",
								 snapshot_check.desired_value, macrostate["total_strain"]));

		}


		/* ---- dyanmics ----- */
		dynamics::extremal_dynamics_step(system);
		history.add_macro(system);

		dynamics::relaxation(system, p.sim.fracture_limit, continue_simulation);
		history.add_macro(system);


		continue_simulation(system.macrostate["total_strain"] < p.sim.monitor_limit,
							"total strain limit reached");
	}

	if(p.out.verbosity and omp_get_thread_num() == 0)
		std::cout << continue_simulation << std::endl;

	timer->leave_subsection("Dynamics");


	/////////////////////////////////////
	//////  post-simultaion step ///////
	///////////////////////////////////

	timer->enter_subsection("Post-simulation step");

	std::string filename = p.out.filename + "_" + write::make_filename(p) + ".h5";
	H5::H5File file(p.out.path + "/" + filename, H5F_ACC_TRUNC);
	write::file_attrs(file, p);
	write::evolution_history(file, history);
	write::snapshots(file, "/snapshots", stress_snapshots);
	file.close();

	/* ---- delete dynamically-allocated objects ----- */
	for(auto &element : elements)
		delete element;

	timer->leave_subsection("Post-simulation step");

	if(p.out.verbosity and omp_get_thread_num() == 0)
		timer->print_summary();

} // run

}; // simulation launcher


} // example



int main(int argc, char *argv[])
{
	cli::Parser parser(argc, argv);
	parser.set_optional<std::string>("f", "file", "./default.cfg",
									 "Name of the input configuration file");
	parser.run_and_exit_if_error();

	dealii::deallog.depth_console(0);


	example::parameters::Parameters p;

	try
	{
		p.load_file(parser.get<std::string>("f"));
	}
	catch(dealii::PathSearch::ExcFileNotFound &)
	{
		p.generate_file(parser.get<std::string>("f"));
		std::cout << "Configuration file " << parser.get<std::string>("f") << " created" << std::endl;
	}

	example::Simulation<2> sim;
	int result = sim.run(p);

	return result;
}