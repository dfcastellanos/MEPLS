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

#include "espci.h"

#if defined(OPENMP)
#include <omp.h>

#endif

#include <cmdparser.hpp>
#include <random>


namespace espci
{

using namespace mepls;


void run(const parameters::Parameters &p, dealii::ConditionalOStream & cout)
{
	// TODO here, dim=2 at least for how we set the average pressure

	/////////////////////////////
	//////  setup system ///////
	///////////////////////////

	constexpr unsigned int dim = 2;

	utils::ContinueSimulation continue_simulation;

	auto timer = utils::TimerSingleton::getInstance();
	timer->enter_subsection("Setting up");

	std::mt19937 generator(p.sim.seed);

	std::vector<espci::element::Anisotropic<dim> *> elements_espci = create_elements<dim>(p,
																						  generator);
	mepls::element::Vector<dim> elements;
	for(auto &element : elements_espci)
		elements.push_back(element);

	snapshot::Check<dim> snapshot_check(p.out.snapshots_min, p.out.snapshots_max,
										p.out.snapshots_interval, p.out.snapshots_sensitivity);
	std::vector<snapshot::Threshold<dim> > threshold_snapshots;
	std::vector<snapshot::Stress<dim> > stress_snapshots;
	std::vector<snapshot::DefGrad<dim> > def_grad_snapshots;
	std::vector<patches::PatchPropertiesSnapshot<dim> > patch_prop_snapshots;
	std::vector<GlobalPropertiesSnapshot<dim>> global_properties_snapshots;
	auto patch_to_element_map = patches::make_patch_to_element_map(elements, p.sim.N_patch_list,
																   p.sim.Nx, p.sim.Ny);
	bool parent_liquid = p.sim.parent_liquid;
	bool thermal_relaxation = p.sim.thermal_relaxation;
	bool reload = p.sim.reload;
	bool het_elasticity = p.sim.het_elasticity;
	bool precalculate = not het_elasticity; // precalculate stress factors when doing patch tests;
	// only valid for homogeneous elasticity
	bool do_ee = p.sim.do_ee; // relax from oi state to ee when doing patch tests

	std::vector<double> theta_list;
	if(p.sim.n_theta == 4)
	{
		theta_list.push_back(0 * M_PI / 180.);
		theta_list.push_back(90 * M_PI / 180.);
		theta_list.push_back(50 * M_PI / 180.);
		theta_list.push_back(140 * M_PI / 180.);
	}
	else if(p.sim.n_theta == 2)
	{
		theta_list.push_back(0 * M_PI / 180.);
		theta_list.push_back(90 * M_PI / 180.);
	}
	else
	{
		for(double z = 0.; z < M_PI; z += M_PI / double(p.sim.n_theta))
			theta_list.push_back(z);
	}

	if(not het_elasticity)
	{
		// If the properties are homogeneous, set them here since the elements,
		// by default, get their properties from a pdf upon their construction.
		// (even if the disorder is very small, this procedure ensures a perfect
		// homogeneity)
		auto CC = utils::tensor::make_mandel_anisotropic_stiffness<dim>(p.mat.average_K_quench, p.mat.average_G_quench,
																		p.mat.average_G_quench, 0.);
		for(auto &element : elements)
			element->C(utils::tensor::mandel_to_standard_rank4<dim>(CC));
	}

	elasticity_solver::LeesEdwards<dim> solver(p.sim.Nx, p.sim.Ny);
	M_Assert(elements.size() == solver.get_n_elements(), "Numbers of mesoscale and finite "
													  "elements do not match");

	// initial (default) assembly using the quench elastic properties.
	for(auto &element : elements)
		solver.set_elastic_properties(element->number(), element->C());
	solver.setup_and_assembly();

	if(het_elasticity)
		mepls::element::calculate_local_stress_coefficients(elements, solver);
	else
		mepls::element::calculate_local_stress_coefficients_central(elements, solver);

	mepls::element::calculate_ext_stress_coefficients(elements, solver);

	system::Standard<dim> system(elements, solver, generator);
	MacroState<dim> &macrostate = system.macrostate;

	if(p.mat.init_eigenstrain)
		apply_initial_eigenstrain<dim>(system, p);

	timer->leave_subsection("Setting up");


	///////////////////////////////////////////
	///////  simulate parent liquid //////////
	/////////////////////////////////////////

	timer->enter_subsection("Simulate parent liquid");

	history::History<dim> liquid_history("parent_liquid");
	system.set_history(liquid_history);
	if(parent_liquid)
	{
		for(auto &element : elements_espci)
		{
			auto conf = element->config();
			conf.temperature = p.mat.temperature_liquid;
			element->config(conf);
		}

		simulate_parent_liquid_KMC(system, liquid_history, p, continue_simulation);
//		simulate_parent_liquid_MH(system, liquid_history, p, continue_simulation);

		// apply instantaneous quench
		for(auto &element : elements_espci)
		{
			auto conf = element->config();
			conf.temperature = 0;
			element->config(conf);
		}
		dynamics::relaxation(system, continue_simulation);

		liquid_history.add_macro(system);
		system.macrostate.clear();

	}

	timer->leave_subsection("Simulate parent liquid");


	////////////////////////////
	//////      AQS     ///////
	//////////////////////////

	history::History<dim> aqs_history("AQS");
	system.set_history(aqs_history);

	/* ----- snapshots ----- */
	timer->enter_subsection("Taking snapshots");

	if(p.out.snapshots.find("slip_thresholds") != std::string::npos)
		threshold_snapshots.push_back(
			snapshot::Threshold(system, p.sim.monitor_name, 0.,
								macrostate[p.sim.monitor_name]));
	if(p.out.snapshots.find("stress") != std::string::npos)
		stress_snapshots.push_back(
			snapshot::Stress(system, p.sim.monitor_name, 0.,
							 macrostate[p.sim.monitor_name]));
	if(p.out.snapshots.find("def_grad") != std::string::npos)
		def_grad_snapshots.push_back(
			snapshot::DefGrad(system, p.sim.monitor_name, 0.,
							  macrostate[p.sim.monitor_name]));

	if(p.out.snapshots.find("patches") != std::string::npos)
		for(auto n_patch : p.sim.N_patch_list)
		{
			cout << ">>> " << n_patch << std::endl;
			patch_prop_snapshots.push_back(
				patches::PatchPropertiesSnapshot<dim>(system,
													  p.sim.monitor_name, 0.,
													  macrostate[p.sim.monitor_name],
													  n_patch, theta_list, precalculate,
													  do_ee));
		}

	if(p.out.snapshots.find("global_properties") != std::string::npos)
		global_properties_snapshots.push_back(
			GlobalPropertiesSnapshot<dim>(system, solver, aqs_history.index(), p.sim.monitor_name,
										  0., macrostate[p.sim.monitor_name]));
	timer->leave_subsection("Taking snapshots");


	/* ----- save struct. properties ----- */
	std::vector<event::RenewSlip<dim>> renewal_vector;
	for(auto &element : elements)
		element->record_structural_properties(renewal_vector);
	aqs_history.add(renewal_vector);


	/* ---- simulation loop ----- */
	timer->enter_subsection("Running AQS");

	continue_simulation(system.macrostate[p.sim.monitor_name] < p.sim.monitor_limit,
						p.sim.monitor_name + " limit reached");


	// adding an (empty) driving event will record the average prestress
	// (which can be slightly non-zero due to numerical accuray) as the
	// first event in the history
	event::Driving<dim> prestress_event;
	prestress_event.activation_protocol = mepls::dynamics::Protocol::prestress;
	system.add(prestress_event);

	aqs_history.add_macro(system);

	double G_old = 0.;
	while(continue_simulation())
	{
		// reassemble the elastic properties if the change if the shear
		// modulus is big enough, if it's different enough from the stationary
		// value and only if we are using homogeneous elastic properties
		double x_global = std::exp(-macrostate["av_vm_plastic_strain"]/p.mat.gamma_pl_trans);
		double G_stat = p.mat.average_G;
		double G_quench = p.mat.average_G_quench;
		double G_new = (G_quench-G_stat)*x_global + G_stat;
		if( not het_elasticity and std::abs(G_new/G_old-1)>0.001 and std::abs(G_old/G_stat-1)>0.005 )
		{
			timer->leave_subsection("Running AQS");
			timer->enter_subsection("Reassembling with new stiffness");

			G_old = G_new;

			double K_stat = p.mat.average_K;
			double K_quench = p.mat.average_K_quench;
			double K_new = (K_quench-K_stat)*x_global + K_stat;

			auto CC = utils::tensor::make_mandel_anisotropic_stiffness<dim>(K_new, G_new,
																			G_new, 0.);
			auto C = utils::tensor::mandel_to_standard_rank4<dim>(CC);

			for(auto &element : elements)
				element->C(C);

			solver.reassemble_with_new_stiffness(C);

			// we call add which will inform the macrostate and the elements
			// about the change in external stress due to the change in global
			// stiffness and will also record that change in the history
			const std::vector<event::Plastic<dim>> added_yielding;
			event::Driving<dim> driving_event_variation_stiffness;
			driving_event_variation_stiffness.activation_protocol = mepls::dynamics::Protocol::variation_stiffness;
			system.add(driving_event_variation_stiffness);

			auto state = solver.get_state();
			mepls::element::calculate_local_stress_coefficients_central(elements, solver);
			mepls::element::calculate_ext_stress_coefficients(elements, solver);
			solver.set_state(state);

			// relax to ensure there are no unstable elements after the change
			// in the elastic properties (if the shear modulus rises with strain,
			// that will lead to stress rises that can unstabilise elements)
			dynamics::relaxation(system, continue_simulation);

			timer->leave_subsection("Reassembling with new stiffness");
			timer->enter_subsection("Running AQS");
		}

		cout << aqs_history.index() << " | " << std::fixed << macrostate["total_strain"]
				  << " " << macrostate["ext_stress"] << std::endl;

		if(snapshot_check(macrostate[p.sim.monitor_name]))
		{
			timer->leave_subsection("Running AQS");
			timer->enter_subsection("Taking snapshots");

			if(p.out.snapshots.find("slip_thresholds") != std::string::npos)
				threshold_snapshots.push_back(
					snapshot::Threshold(system, p.sim.monitor_name,
										snapshot_check.desired_value,
										macrostate[p.sim.monitor_name]));
			if(p.out.snapshots.find("stress") != std::string::npos)
				stress_snapshots.push_back(
					snapshot::Stress(system, p.sim.monitor_name,
									 snapshot_check.desired_value, macrostate[p.sim.monitor_name]));
			if(p.out.snapshots.find("def_grad") != std::string::npos)
				def_grad_snapshots.push_back(
					snapshot::DefGrad(system, p.sim.monitor_name,
									  snapshot_check.desired_value,
									  macrostate[p.sim.monitor_name]));
			if(p.out.snapshots.find("patches") != std::string::npos)
				for(auto n_patch : p.sim.N_patch_list)
				{
					cout << ">>> " << n_patch << std::endl;
					patch_prop_snapshots.push_back(
						patches::PatchPropertiesSnapshot<dim>(system,
															  p.sim.monitor_name,
															  snapshot_check.desired_value,
															  macrostate[p.sim.monitor_name],
															  n_patch,
															  theta_list, precalculate, do_ee));
				}

			timer->leave_subsection("Taking snapshots");
			timer->enter_subsection("Running AQS");
		}


		/* ---- dyanmics ----- */
		dynamics::finite_extremal_dynamics_step(1e-4 * 0.5, system);
		aqs_history.add_macro(system);

		dynamics::relaxation(system, continue_simulation);
		aqs_history.add_macro(system);


		continue_simulation(system.macrostate[p.sim.monitor_name] < p.sim.monitor_limit,
							p.sim.monitor_name + " limit reached");
	}

	auto solver_state_end_AQS = solver.get_state();
	auto macrostate_end_AQS = system.macrostate;

	cout << continue_simulation << std::endl;

	timer->leave_subsection("Running AQS");


	//////////////////////////////////
	//////  thermal relaxation //////
	////////////////////////////////

	timer->enter_subsection("Thermal relaxation");

	history::History<dim> thermal_relaxation_hist("thermal_relaxation");

	if(thermal_relaxation)
	{
		mepls::element::Vector<dim> elements_replica;
		for(auto &element : elements_espci)
		{
			auto element_copy = element->make_copy();

			auto conf = element_copy->config();
			conf.temperature = p.mat.temperature_relaxation;
			element_copy->config(conf);

			elements_replica.push_back( element_copy );
		}

		// the solver state is copied, so it contains the eigenstrain and the load
		// in this case the elements elastic field need not be converted into a
		// prestress before
		solver.set_state(solver_state_end_AQS);

		auto system_replica = system.get_new_instance(elements_replica, solver,
													  system.generator);
		system_replica->macrostate = macrostate_end_AQS;
		system_replica->set_history(thermal_relaxation_hist);

		// initiate dynamics using copied system
		mepls::dynamics::KMC<dim> kmc;
		mepls::utils::ContinueSimulation continue_relaxing;
		auto &macrostate = system_replica->macrostate;
		thermal_relaxation_hist.add_macro( *system_replica );

		while(continue_relaxing())
		{

			cout << thermal_relaxation_hist.index() << " | " << std::fixed
					  << macrostate["total_strain"]
					  << " " << macrostate["ext_stress"]
					  << " " << macrostate["pressure"]
					  << " " << macrostate["time"]
					  << std::endl;

			kmc(*system_replica);
			thermal_relaxation_hist.add_macro( *system_replica );

			mepls::dynamics::relaxation(*system_replica, continue_relaxing);
			thermal_relaxation_hist.add_macro( *system_replica );

			continue_relaxing(macrostate["ext_stress"] > 0, "System relaxed");
		}

		delete system_replica;

		cout << continue_relaxing << std::endl;
		cout << "Relaxing finished" << std::endl;
	}

	timer->leave_subsection("Thermal relaxation");

	///////////////////////////////
	//////  reloading step ///////
	/////////////////////////////

	history::History<dim> aqs_unloading("AQS_unloading");
	history::History<dim> aqs_reload_forward("AQS_reload_forward");
	history::History<dim> aqs_reload_backward("AQS_reload_backward");

	if(reload)
	{
		timer->enter_subsection("Reloading");

		// the solver state is copied, so it contains the eigenstrain and the load
		// in this case the elements elastic field need not be converted into a
		// prestress before
		solver.set_state(solver_state_end_AQS);

		system.macrostate = macrostate_end_AQS;

		system.set_history(aqs_unloading);

		utils::ContinueSimulation continue_unloading;
		auto &macrosate = system.macrostate;
		aqs_unloading.add_macro(system);

		while(continue_unloading())
		{
			cout << aqs_unloading.index() << " | " << std::fixed
					  << macrostate["total_strain"] << " " << macrostate["ext_stress"] << " "
					  << macrostate["pressure"] << std::endl;

			dynamics::finite_extremal_dynamics_step(1e-4 * 0.5, system, false);
			aqs_unloading.add_macro(system);

			dynamics::relaxation(system, continue_unloading);
			aqs_unloading.add_macro(system);

			continue_unloading(macrostate["ext_stress"] > 0, "System unloaded");
		}

		cout << continue_unloading << std::endl;

		perform_reloading(system, aqs_reload_forward, true, p);
		perform_reloading(system, aqs_reload_backward, false, p);

		timer->leave_subsection("Reloading");
	}


	/////////////////////////////////////
	//////  post-simultaion step ///////
	///////////////////////////////////

	std::string filename = p.out.filename + "_" + write::make_filename(p) + ".h5";
	H5::H5File file(p.out.path + "/" + filename, H5F_ACC_TRUNC);
	write::file_attrs(file, p);

	if(parent_liquid)
		write::evolution_history(file, liquid_history);
	if(thermal_relaxation)
		write::evolution_history(file, thermal_relaxation_hist);
	write::evolution_history(file, aqs_history);
	if(reload)
	{
		write::evolution_history(file, aqs_unloading);
		write::evolution_history(file, aqs_reload_forward);
		write::evolution_history(file, aqs_reload_backward);
	}

	write::patch_info<dim>(file, "/patch_info", patch_to_element_map);
	write::element_info(file, "/element_info", elements);
	write::snapshots(file, "/snapshots", threshold_snapshots, stress_snapshots, def_grad_snapshots,
					 patch_prop_snapshots, global_properties_snapshots);
	file.close();


	/* ---- delete dynamically-allocated objects ----- */
	for(auto &element : elements)
		delete element;

	if(p.out.verbosity and omp_get_thread_num() == 0)
		timer->print_summary();
} // run


} // espci


int main(int argc, char *argv[])
{
	cli::Parser parser(argc, argv);
	parser.set_optional<std::string>("f", "file", "./default.cfg",
									 "Name of the input configuration file");
	parser.run_and_exit_if_error();

	dealii::deallog.depth_console(0);


	espci::parameters::Parameters p;

	try
	{
		p.load_file(parser.get<std::string>("f"));
	}
	catch(dealii::PathSearch::ExcFileNotFound &)
	{
		p.generate_file(parser.get<std::string>("f"));
		std::cout << "Configuration file " << parser.get<std::string>("f") << " created" << std::endl;
	}


	unsigned int n_rep = p.sim.n_rep;
	if(n_rep < omp_get_max_threads())
		n_rep = omp_get_max_threads();

	// initialize the master engine with the master seed
	std::srand(p.sim.seed);

	#pragma omp parallel
	{
		unsigned int n_threads = omp_get_max_threads();
		unsigned int id = omp_get_thread_num();
		unsigned int rep_per_thread = int( n_rep / n_threads );

		dealii::ConditionalOStream cout(std::cout, id==0 and p.out.verbosity);

		for(unsigned int n = 0; n < rep_per_thread; ++n)
		{
			espci::parameters::Parameters p_thread = p;
			#pragma critical
			{
				// generate a seed for a simulation run
				p_thread.sim.seed = std::rand();
			};

			espci::run(p_thread, cout);
		}

	}

}