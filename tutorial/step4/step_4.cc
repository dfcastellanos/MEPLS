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
#include <mepls/history.h>

// MEPLS built-in dynamics
#include <mepls/dynamics.h>

// to parse command line arguments
#include <cmdparser.hpp>

// to parse input parameters files
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/path_search.h>

// to save output data in JSON format
#include <boost/property_tree/json_parser.hpp>


struct Parameters
{
	// these are the parameters and we will use in the simulation, initialized with some default values
	unsigned int seed = 1234567;
	unsigned int Nx = 32;
	unsigned int Ny = 32;
	double G = 30.;
	double nu = 0.3;
	double gamma = 0.05;
	double k = 6.;
	double strain_limit = 0.05;
	double lambda_init = 1.;
	double lambda_renew = 1.;
	std::string filename = "out.json";

	void declare_entries(dealii::ParameterHandler &prm)
	{
		// We declare the entries of the parameters text file. Each entry matches the name of a
		// simulation parameters (although it doesn't need to) and has a default value. We use as
		// the default value is the same value of the variables declared above

		prm.enter_subsection("Section1");

		prm.declare_entry("seed", mepls::utils::str::to_string(seed), dealii::Patterns::Integer(0),
						  "");
		prm.declare_entry("Nx", mepls::utils::str::to_string(Nx), dealii::Patterns::Integer(0), "");
		prm.declare_entry("Ny", mepls::utils::str::to_string(Ny), dealii::Patterns::Integer(0), "");
		prm.declare_entry("G", mepls::utils::str::to_string(G), dealii::Patterns::Double(0.0), "");
		prm.declare_entry("nu", mepls::utils::str::to_string(nu), dealii::Patterns::Double(0.0),
						  "");
		prm.declare_entry("gamma", mepls::utils::str::to_string(gamma),
						  dealii::Patterns::Double(0.0), "");
		prm.declare_entry("strain_limit", mepls::utils::str::to_string(strain_limit),
						  dealii::Patterns::Double(0.0), "");
		prm.declare_entry("lambda_init", mepls::utils::str::to_string(lambda_init),
						  dealii::Patterns::Double(0.0), "");
		prm.declare_entry("lambda_renew", mepls::utils::str::to_string(lambda_renew),
						  dealii::Patterns::Double(0.0), "");
		prm.declare_entry("k", mepls::utils::str::to_string(k), dealii::Patterns::Double(0.0), "");
		prm.declare_entry("filename", filename, dealii::Patterns::FileName(), "");

		prm.leave_subsection();
	}

	void load_entries(dealii::ParameterHandler &prm)
	{
		prm.enter_subsection("Section1");

		// We define how each parameter gets its value from a parameters file entry

		seed = prm.get_integer("seed");
		Nx = prm.get_integer("Nx");
		Ny = prm.get_integer("Ny");
		G = prm.get_double("G");
		nu = prm.get_double("nu");
		gamma = prm.get_double("gamma");
		k = prm.get_double("k");
		strain_limit = prm.get_double("strain_limit");
		lambda_init = prm.get_double("lambda_init");
		lambda_renew = prm.get_double("lambda_renew");
		filename = prm.get("filename");

		prm.leave_subsection();
	}

	void load_file(const std::string &filename)
	{
		// load a parameters file

		dealii::ParameterHandler prm;
		declare_entries(prm);
		prm.parse_input(filename);
		load_entries(prm);
	}

	void generate_file(const std::string &filename)
	{
		// generate a parameters file template

		std::ofstream outfile(filename);
		dealii::ParameterHandler prm;
		declare_entries(prm);
		prm.print_parameters(outfile, dealii::ParameterHandler::OutputStyle::Text);
	}

};


template<int dim>
void write_data(
	const mepls::history::History<dim> &sim_history, const Parameters &p)
{
	boost::property_tree::ptree data_tree;

	// write some metadata, such as a simulation name and description
	data_tree.put("Name", "Step4");
	data_tree.put("Description", "System driven in the athermal quasistatic limit");

	// write the simulation parameters
	data_tree.put("Parameters.dim", 2);
	data_tree.put("Parameters.seed", p.seed);
	data_tree.put("Parameters.Nx", p.Nx);
	data_tree.put("Parameters.Ny", p.Ny);
	data_tree.put("Parameters.G", p.G);
	data_tree.put("Parameters.nu", p.nu);
	data_tree.put("Parameters.gamma", p.gamma);
	data_tree.put("Parameters.lambda_renew", p.lambda_renew);
	data_tree.put("Parameters.lambda_init", p.lambda_init);
	data_tree.put("Parameters.k", p.k);
	data_tree.put("Parameters.strain_limit", p.strain_limit);


	// write some metadata, such as a simulation name and description
	data_tree.put("Name", "Step4");
	data_tree.put("Description", "System driven in the athermal quasistatic limit");

	// write the simulation parameters
	data_tree.put("Parameters.dim", 2);
	data_tree.put("Parameters.seed", p.seed);
	data_tree.put("Parameters.Nx", p.Nx);
	data_tree.put("Parameters.Ny", p.Ny);
	data_tree.put("Parameters.G", p.G);
	data_tree.put("Parameters.nu", p.nu);
	data_tree.put("Parameters.gamma", p.gamma);
	data_tree.put("Parameters.lambda_renew", p.lambda_renew);
	data_tree.put("Parameters.lambda_init", p.lambda_init);
	data_tree.put("Parameters.k", p.k);
	data_tree.put("Parameters.strain_limit", p.strain_limit);

	// We write the event histories with CSV format to a string.
	// Here, we write only the columns of interest, but there are more available (see the
	// documentation).
	std::ostringstream plastic_events_csv;
	plastic_events_csv << "index,element,eigenstrain_00,eigenstrain_11,eigenstrain_01\n";
	for(auto &row : sim_history.plastic)
		plastic_events_csv << row.index << "," << row.element << "," << row.eigenstrain_00 << ","
						   << row.eigenstrain_11 << "," << row.eigenstrain_01 << "\n";

	data_tree.put("Data.plastic_events", plastic_events_csv.str());


	std::ostringstream driving_events_csv;
	driving_events_csv << "index,dext_stress,dtotal_strain\n";
	for(auto &row : sim_history.driving)
		driving_events_csv << row.index << "," << row.dext_stress << "," << row.dtotal_strain
						   << "\n";

	data_tree.put("Data.driving_events", driving_events_csv.str());


	std::ostringstream macro_evolution_csv;
	macro_evolution_csv << "index,ext_stress,total_strain,time,av_vm_stress,av_vm_plastic_strain\n";
	for(auto &row : sim_history.macro_evolution)
		macro_evolution_csv << row.index << "," << row.ext_stress << "," << row.total_strain << ","
							<< row.time << "," << row.av_vm_stress << ","
							<< row.av_vm_plastic_strain << "\n";

	data_tree.put("Data.macro_evolution", macro_evolution_csv.str());


	// write the data tree to a JSON file
	std::ofstream output_file(p.filename);
	boost::property_tree::json_parser::write_json(output_file, data_tree);
	output_file.close();
}


void run(const Parameters &p)
{
	//----- SETUP ------

	// we do the same as in the previous tutorial, but this time we use the
	// parameters from the input Parameters object

	constexpr unsigned dim = 2;
	std::mt19937 generator(p.seed);

	dealii::SymmetricTensor<4, dim> C = mepls::utils::tensor::make_isotropic_stiffness<dim>(p.G,
																							p.nu);

	mepls::element::Vector<dim> elements;

	// we also store the elements into a vector that knows the derived class, so we can access
	// the example::element::Scalar<dim>::conf struct later (this member cannot be accessed
	// through the pointers to the base mepls::element::Element<dim>)
	std::vector<example::element::Scalar<dim> *> elements_scalar;

	for(double n = 0; n < p.Nx * p.Ny; ++n)
	{
		example::element::Scalar<dim>::Config conf;
		conf.number = n;
		conf.gamma = p.gamma;
		conf.lambda = p.lambda_init;
		conf.k = p.k;

		auto element = new example::element::Scalar<dim>(conf, generator);
		element->C(C);

		elements_scalar.push_back(element);
		elements.push_back(element);
	}

	mepls::elasticity_solver::LeesEdwards<dim> solver(p.Nx, p.Ny);
	for(auto &element : elements)
		solver.set_elastic_properties(element->number(), element->C());
	solver.setup_and_assembly();

	mepls::element::calculate_ext_stress_coefficients(elements, solver);
	mepls::element::calculate_local_stress_coefficients_central(elements, solver);

	mepls::system::Standard<dim> system(elements, solver, generator);

	mepls::history::History<dim> sim_history("Simulation_history");

	system.set_history(sim_history);

	// when the elements were created, their slip systems were initialized with thresholds from a
	// Weibull distribution with scale lambda_init. We change it now to lambda_renew. Therefore,
	// when plastic deformation occurs, the renew local thresholds will have a different
	// average
	for(auto &element : elements_scalar)
		element->conf.lambda = p.lambda_renew;


	//----- DYNAMICS ------

	// This object will allow us to check for different conditions by which
	// the simulation might stop. Different conditions might be checked, but as long as one of
	// them evaluates to false, when continue_simulation() is called it will return false
	mepls::utils::ContinueSimulation continue_simulation;

	sim_history.add_macro(system);

	// run while the continue_simulation object says so
	while(continue_simulation())
	{
		std::cout << system.macrostate["total_strain"] << " " << system.macrostate["ext_stress"]
				  << std::endl;


		// these are the same dynamics we implemented in the previous tutorial, but using MEPLS built-in

		// apply an external strain increment of 0.01%
		mepls::dynamics::finite_extremal_dynamics_step(1e-4, system);
		sim_history.add_macro(system);

		// perform and avalanche of slip events. By passing the continue_simulation object,
		// the relaxation function can set its own condition for stopping the simulation.
		// Specifically, it will check if the avalanche size overcomes a certain maximum upper limit
		mepls::dynamics::relaxation(system, continue_simulation);
		sim_history.add_macro(system);


		// check if the strain has reached the strain limit. If it has, the next time continue_simulation() is called
		// it will return false, so we will exit the main loop
		continue_simulation(system.macrostate["total_strain"] < p.strain_limit,
							"total strain limit reached");

	}

	// print the message of the stopping condition that was met
	std::cout << continue_simulation << std::endl;

	for(auto &element : elements)
		delete element;


	write_data(sim_history, p);
}


int main(int argc, char *argv[])
{
	// Read the command line arguments. We define the -f 'filename' to pass the
	// path to the parameters file
	cli::Parser parser(argc, argv);
	parser.set_optional<std::string>("f", "file", "./default.prm",
									 "Name of the input configuration file");
	parser.run_and_exit_if_error();
	// you can check https://github.com/FlorianRappl/CmdParser for a further documentation of cli::Parser



	// Create the parameters object
	Parameters p;

	// We try to load the parameters file, but if it doesn't exist, we generate a new one with the name
	// default.prm and default values
	try
	{
		p.load_file(parser.get<std::string>("f"));
	}
	catch(dealii::PathSearch::ExcFileNotFound &)
	{
		p.generate_file(parser.get<std::string>("f"));
		std::cout << "Configuration file " << parser.get<std::string>("f") << " created"
				  << std::endl;
		return 1;
	}


	// run the simulation
	run(p);


	return 0;
}
