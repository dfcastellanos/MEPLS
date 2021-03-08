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

template<int dim>
void run(const parameters::Standard &p)
{

	// TODO here, dim=2 at least for how we set the average pressure

	/////////////////////////////
	//////  setup system ///////
	///////////////////////////

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
	auto patch_to_element_map = patches::make_patch_to_element_map(elements, p.sim.N_probe_list,
																   p.sim.Nx, p.sim.Ny);
	bool kmc_quench = p.sim.kmc_quench;
	bool kmc_relaxation = p.sim.kmc_relaxation;
	bool reload = p.sim.reload;
	bool het_elasticity = p.sim.het_elasticity;
	bool precalculate = not het_elasticity; // precalculate stress factors when doing local probes; only valid for homogeneous elasticity
	bool do_ee = p.sim.do_ee; // relax from oi state to ee when doing local probes

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
		// make sure that properties are heterogeneous, even if we set weibull shapes to very big values,
		// this ensures a perfect homogeneity
		auto CC = utils::tensor::make_mandel_anisotropic_stiffness<dim>(p.mat.average_K, p.mat.average_G,
																		p.mat.average_G, 0.);
		for(auto &element : elements)
			element->C(utils::tensor::mandel_to_standard_rank4<dim>(CC));
	}

	// assemble solver to be used in the simulation
	elasticity_solver::LeesEdwards<dim> solver(p.sim.Nx, p.sim.Ny);
	M_Assert(elements.size() == solver.get_n_elements(), "Numbers of mesoscale and finite "
													  "elements do not match");
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

	timer->leave_subsection("Setting up");


	/////////////////////////////////
	///////  quench state //////////
	///////////////////////////////

	timer->enter_subsection("Creating quench");

	std::normal_distribution<double> pressure_average_dist(0., p.mat.prestress_std_av_pressure);
	double av_pressure_fluctuation = pressure_average_dist(generator);
	for(auto &element : elements_espci)
	{
		//		dealii::SymmetricTensor<2,dim> av_pressure;
		//		av_pressure[0][0] = p.mat.prestress_av_pressure + av_pressure_fluctuation;
		//		av_pressure[1][1] = p.mat.prestress_av_pressure + av_pressure_fluctuation;
		//		element->prestress( element->prestress() + av_pressure);

			auto conf = element->config();
			conf.temperature = p.mat.temperature_liquid;
			element->config(conf);
	}


	//   quench::equilibrate_initial_structure_rejection(system, p, generator);

	history::EventAndMacro<dim> kmc_history("KMC_quench");
	system.set_history(kmc_history);
	if(kmc_quench)
		quench::run_thermal_evolution(system, kmc_history, p, continue_simulation);
	quench::convert_state_to_quench(system);

	//	Set the threshold renewal properties before performing any quench-state shear tests
	auto weibull = [&](double x, double k, double lambda)
	{ return std::pow(x, k - 1) * std::exp(-std::pow(x / lambda, k)); };
	auto func = std::bind(weibull, std::placeholders::_1, p.mat.k, p.mat.lambda);
	auto threshold_distribution_ptr = utils::rand::create_distribution(func, 1e-3, 0., 1e-7, 8);

	for(auto &element : elements_espci)
	{
		auto conf = element->config();
		conf.temperature = 0;
		element->config(conf);

		element->threshold_distribution(
			threshold_distribution_ptr); // next call to renew_stuctrua_properties will use this distribution
	}

	timer->leave_subsection("Creating quench");


	////////////////////////////
	//////      AQS     ///////
	//////////////////////////

	history::EventAndMacro<dim> aqs_history("AQS");
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

	if(p.out.snapshots.find("local_probe") != std::string::npos)
		for(auto n_probe : p.sim.N_probe_list)
		{
			if(p.out.verbosity and omp_get_thread_num() == 0)
				std::cout << ">>> " << n_probe << std::endl;
			patch_prop_snapshots.push_back(
				patches::PatchPropertiesSnapshot<dim>(system,
													  p.sim.monitor_name, 0.,
													  macrostate[p.sim.monitor_name],
													  n_probe, theta_list, precalculate,
													  do_ee));
		}

	if(p.out.snapshots.find("global_properties") != std::string::npos)
		global_properties_snapshots.push_back(
			GlobalPropertiesSnapshot<dim>(system, solver, aqs_history.index, p.sim.monitor_name,
										  0., macrostate[p.sim.monitor_name]));
	timer->leave_subsection("Taking snapshots");


	/* ----- save struct. properties ----- */
	std::vector<event::RenewSlip<dim>> renewal_vector;
	for(auto &element : elements)
		element->record_structural_properties(renewal_vector);
	aqs_history.add(renewal_vector);


	/* ---- simulation loop ----- */
	timer->enter_subsection("Running AQS");

	dynamics::relaxation(system, p.sim.fracture_limit, continue_simulation);


	continue_simulation(system.macrostate[p.sim.monitor_name] < p.sim.monitor_limit,
						p.sim.monitor_name + " limit reached");

	while(continue_simulation())
	{
		if(p.out.verbosity and omp_get_thread_num() == 0)
			std::cout << aqs_history.index << " | " << std::fixed << macrostate["total_strain"]
					  << " " << macrostate["ext_stress"] << " " << macrostate["pressure"]
					  << std::endl;

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
			if(p.out.snapshots.find("local_probe") != std::string::npos)
				for(auto n_probe : p.sim.N_probe_list)
				{
					if(p.out.verbosity and omp_get_thread_num() == 0)
						std::cout << ">>> " << n_probe << std::endl;
					patch_prop_snapshots.push_back(
						patches::PatchPropertiesSnapshot<dim>(system,
															  p.sim.monitor_name,
															  snapshot_check.desired_value,
															  macrostate[p.sim.monitor_name],
															  n_probe,
															  theta_list, precalculate, do_ee));
				}

			timer->leave_subsection("Taking snapshots");
			timer->enter_subsection("Running AQS");
		}

		aqs_history.add_macro(system);

		/* ---- dyanmics ----- */
		dynamics::finite_extremal_dynamics_step(1e-4 * 0.5, system);
		dynamics::relaxation(system, p.sim.fracture_limit, continue_simulation);

		continue_simulation(system.macrostate[p.sim.monitor_name] < p.sim.monitor_limit,
							p.sim.monitor_name + " limit reached");
	}

	auto solver_state_end_AQS = solver.get_state();
	auto macrostate_end_AQS = system.macrostate;

	if(p.out.verbosity and omp_get_thread_num() == 0)
		std::cout << continue_simulation << std::endl;

	timer->leave_subsection("Running AQS");


	//////////////////////////////////
	//////  thermal relaxation //////
	////////////////////////////////

	timer->enter_subsection("Thermal relaxation");

	history::EventAndMacro<dim> kmc_relaxation_hist("KMC_relaxation");

	if(kmc_relaxation)
	{
		mepls::element::Vector<dim> elements_replica;
		for(auto &element : elements_espci)
			elements_replica.push_back( element->make_copy() );

		// the solver state is copied, so it contains the eigenstrain and the load
		// in this case the elements elastic field need not be converted into a
		// prestress before
		solver.set_state(solver_state_end_AQS);

		auto system_replica = system.get_new_instance(elements_replica, solver,
													  system.generator);
		system_replica->macrostate = macrostate_end_AQS;
		system_replica->set_history(kmc_relaxation_hist);

		// initiate dynamics using copied system
		mepls::dynamics::KMC<dim> kmc;
		mepls::utils::ContinueSimulation continue_relaxing;
		auto &macrostate = system_replica->macrostate;
		while(continue_relaxing())
		{
			// modify the temp. with time as in the MD system
			double t = macrostate["time"];
			double tp = 7.3*100;
			double T = t < tp ? 0.5*p.mat.temperature_relaxation*(t/tp+1) : p.mat.temperature_relaxation;

			// update the elements temperature
			for(auto &element : elements_replica)
			{
				auto element_espci = static_cast<espci::element::Anisotropic<dim> *>(element);
				auto conf = element_espci->config();
				conf.temperature = T;
				element_espci->config(conf);
			}

			if(p.out.verbosity and omp_get_thread_num() == 0)
				std::cout << kmc_relaxation_hist.index << " | " << std::fixed
						  << macrostate["total_strain"]
						  << " " << macrostate["ext_stress"]
						  << " " << macrostate["pressure"]
						  << " " << macrostate["time"]
						  << std::endl;

			kmc(*system_replica, continue_relaxing);
			mepls::dynamics::relaxation(*system_replica, p.sim.fracture_limit, continue_relaxing);

			kmc_relaxation_hist.add_macro( *system_replica );
			// set the value of the temperature in the macro. evol. dataset
			// TODO add_macro( *system_replica ) should have access to the
			// elements' temperature and compute the average
			kmc_relaxation_hist.macro_evolution.back().temperature = T;

			continue_relaxing(std::abs(macrostate["time"]) < 5000, "System relaxed");
		}

		delete system_replica;

		if(p.out.verbosity and omp_get_thread_num() == 0)
		{
			std::cout << continue_relaxing << std::endl;
			std::cout << "Relaxing finished" << std::endl;
		}
	}

	timer->leave_subsection("Thermal relaxation");

	///////////////////////////////
	//////  reloading step ///////
	/////////////////////////////

	history::EventAndMacro<dim> aqs_unloading("AQS_unloading");
	history::EventAndMacro<dim> aqs_reload_forward("AQS_reload_forward");
	history::EventAndMacro<dim> aqs_reload_backward("AQS_reload_backward");

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
		while(continue_unloading())
		{
			if(p.out.verbosity and omp_get_thread_num() == 0)
				std::cout << aqs_unloading.index << " | " << std::fixed
						  << macrostate["total_strain"] << " " << macrostate["ext_stress"] << " "
						  << macrostate["pressure"] << std::endl;

			aqs_unloading.add_macro(system);

			dynamics::finite_extremal_dynamics_step(1e-4 * 0.5, system, false);
			dynamics::relaxation(system, p.sim.fracture_limit, continue_unloading);

			continue_unloading(macrostate["ext_stress"] > 0, "System unloaded");
		}

		if(p.out.verbosity and omp_get_thread_num() == 0)
			std::cout << continue_unloading << std::endl;

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

	if(kmc_quench)
		write::evolution_history(file, kmc_history);
	if(kmc_relaxation)
		write::evolution_history(file, kmc_relaxation_hist);
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
}

} // espci


int main(int argc, char *argv[])
{
	cli::Parser parser(argc, argv);
	parser.set_optional<std::string>("f", "file", "./default.cfg",
									 "Name of the input configuration file");
	parser.run_and_exit_if_error();

	dealii::deallog.depth_console(0);

	try
	{

		espci::parameters::Standard p;

		try
		{
			p.load_file(parser.get<std::string>("f"));

#if defined(OPENMP)

			if(p.sim.n_rep < omp_get_max_threads())
				p.sim.n_rep = omp_get_max_threads();

			// create a seed to initialize the rnd engine of each thread
			srand(p.sim.seed);
			std::vector<unsigned int> seed_for_thread;
			for(unsigned int n = 0; n < omp_get_max_threads(); ++n)
				seed_for_thread.push_back(rand());

#pragma omp parallel
			{
				auto pp = p;
				unsigned int np = omp_get_num_threads();
				unsigned int id = omp_get_thread_num();
				unsigned int nmin = id * std::floor(p.sim.n_rep / omp_get_max_threads());
				unsigned int nmax = (id + 1) * std::floor(p.sim.n_rep / omp_get_max_threads());
				std::mt19937 generator_of_seeds(seed_for_thread[id]);

				for(unsigned int n = nmin; n < nmax; ++n)
				{
					if(id == 0 and p.out.verbosity)
						std::cout << double(n) / double(nmax) * 100 << "%" << std::endl;

					generator_of_seeds.discard(1000);
					pp.sim.seed = generator_of_seeds();

					espci::run<2>(pp);
				}

			} // parallel region

			#else
			espci::run<2>(p);
			#endif // defined(OPENMP)

		}
		catch(dealii::PathSearch::ExcFileNotFound &)
		{
			p.generate_file(parser.get<std::string>("f"));
			std::cout << "Configuration file " << parser.get<std::string>("f") << " created" << std::endl;
		}

	}
	catch(std::exception &exc)
	{
		std::cerr << std::endl << std::endl << "----------------------------------------------------" << std::endl;
		std::cerr << "Exception on processing: " << std::endl << exc.what() << std::endl << "Aborting!" << std::endl
				  << "----------------------------------------------------" << std::endl;
		return 1;
	}
	catch(...)
	{
		std::cerr << std::endl << std::endl << "----------------------------------------------------" << std::endl;
		std::cerr << "Unknown exception!" << std::endl << "Aborting!" << std::endl
				  << "----------------------------------------------------" << std::endl;
		return 1;
	}

	return 0;
}