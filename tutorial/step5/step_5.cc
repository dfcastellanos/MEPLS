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
#include <mepls/dynamics.h>
#include <cmdparser.hpp>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/path_search.h>
#include <boost/property_tree/json_parser.hpp>

// new headers
#include <omp.h>
#include <deal.II/base/conditional_ostream.h>


struct Parameters
{
   unsigned int seed = 1234567;
   unsigned int n_rep = 1;
   unsigned int Nx = 32;
   unsigned int Ny = 32;
   double G = 30.;
   double nu = 0.3;
   double gamma = 0.05;
   double k = 6.;
   double strain_limit_aqs = 0.05;
   double time_limit_creep = 1e5;
   double ext_stress_creep = 0.;
   double lambda = 1.;
   double temperature = 1.;
   bool verbose = true;

   void declare_entries(dealii::ParameterHandler &prm)
   {
      prm.enter_subsection("Section1");

      prm.declare_entry("seed", mepls::utils::str::to_string(seed), dealii::Patterns::Integer(0), "");
      prm.declare_entry("n_rep", mepls::utils::str::to_string(n_rep), dealii::Patterns::Integer(0), "");
      prm.declare_entry("Nx", mepls::utils::str::to_string(Nx), dealii::Patterns::Integer(0), "");
      prm.declare_entry("Ny", mepls::utils::str::to_string(Ny), dealii::Patterns::Integer(0), "");
      prm.declare_entry("G", mepls::utils::str::to_string(G), dealii::Patterns::Double(0.0), "");
      prm.declare_entry("nu", mepls::utils::str::to_string(nu), dealii::Patterns::Double(0.0), "");
      prm.declare_entry("gamma", mepls::utils::str::to_string(gamma), dealii::Patterns::Double(0.0), "");
      prm.declare_entry("strain_limit_aqs", mepls::utils::str::to_string(strain_limit_aqs), dealii::Patterns::Double(0.0), "");
      prm.declare_entry("time_limit_creep", mepls::utils::str::to_string(time_limit_creep), dealii::Patterns::Double(0.0), "");
      prm.declare_entry("ext_stress_creep", mepls::utils::str::to_string(ext_stress_creep), dealii::Patterns::Double(0.0), "");
      prm.declare_entry("temperature", mepls::utils::str::to_string(temperature), dealii::Patterns::Double(0.0), "");
      prm.declare_entry("lambda", mepls::utils::str::to_string(lambda), dealii::Patterns::Double(0.0), "");
      prm.declare_entry("k", mepls::utils::str::to_string(k), dealii::Patterns::Double(0.0), "");
	  prm.declare_entry("verbose", mepls::utils::str::to_string(verbose), dealii::Patterns::Bool(), "");

      prm.leave_subsection();
   }

   void load_entries(dealii::ParameterHandler &prm)
   {
      prm.enter_subsection("Section1");

      seed = prm.get_integer("seed");
      n_rep = prm.get_integer("n_rep");
      Nx = prm.get_integer("Nx");
      Ny = prm.get_integer("Ny");
      G = prm.get_double("G");
      nu = prm.get_double("nu");
      gamma = prm.get_double("gamma");
      k = prm.get_double("k");
      strain_limit_aqs = prm.get_double("strain_limit_aqs");
      time_limit_creep = prm.get_double("time_limit_creep");
      ext_stress_creep = prm.get_double("ext_stress_creep");
      temperature = prm.get_double("temperature");
      lambda = prm.get_double("lambda");
      verbose = prm.get_bool("verbose");

      prm.leave_subsection();
   }

   void load_file(const std::string &filename)
   {
      dealii::ParameterHandler prm;
      declare_entries(prm);
      prm.parse_input(filename);
      load_entries(prm);
   }

   void generate_file(const std::string &filename)
   {
      std::ofstream outfile(filename);
      dealii::ParameterHandler prm;
      declare_entries(prm);
      prm.print_parameters(outfile, dealii::ParameterHandler::OutputStyle::Text);
   }

};


template<int dim>
void write_data(const mepls::history::History<dim> &creep_history,
                const mepls::history::History<dim> &aqs_history,
                const Parameters &p)
{
   boost::property_tree::ptree data_tree;

   data_tree.put("Name", "Step5");
   data_tree.put("Description", "System undergoing creep deformation, and then driven in "
								"athermal quasistatic shear");

   data_tree.put("Parameters.dim", 2);
   data_tree.put("Parameters.seed", p.seed);
   data_tree.put("Parameters.Nx", p.Nx);
   data_tree.put("Parameters.Ny", p.Ny);
   data_tree.put("Parameters.G", p.G);
   data_tree.put("Parameters.nu", p.nu);
   data_tree.put("Parameters.gamma", p.gamma);
   data_tree.put("Parameters.lambda", p.lambda);
   data_tree.put("Parameters.temperature", p.temperature);
   data_tree.put("Parameters.k", p.k);
   data_tree.put("Parameters.ext_stress_creep", p.ext_stress_creep);
   data_tree.put("Parameters.time_limit_creep", p.time_limit_creep);
   data_tree.put("Parameters.strain_limit_aqs", p.strain_limit_aqs);

   	// -------- creep history ----------
	{
	   std::ostringstream plastic_events_csv;
	   plastic_events_csv << "index,element,eigenstrain_00,eigenstrain_11,eigenstrain_01\n";
	   for(auto &row : creep_history.plastic)
		  plastic_events_csv << row.index << "," << row.element << "," << row.eigenstrain_00 << ","
							 << row.eigenstrain_11 << "," << row.eigenstrain_01 << "\n";

	   data_tree.put("Data.creep.plastic_events", plastic_events_csv.str());


	   std::ostringstream driving_events_csv;
	   driving_events_csv << "index,dtime,dext_stress,dtotal_strain\n";
	   for(auto &row : creep_history.driving)
		  driving_events_csv << row.index << "," << row.dtime << ","
							 << row.dext_stress << "," << row.dtotal_strain << "\n";

	   data_tree.put("Data.creep.driving_events", driving_events_csv.str());


	   std::ostringstream macro_evolution_csv;
	   macro_evolution_csv << "index,ext_stress,total_strain,time,av_vm_stress,av_vm_plastic_strain\n";
	   for(auto &row : creep_history.macro_evolution)
		  macro_evolution_csv << row.index << "," << row.ext_stress << ","
							 << row.total_strain << "," << row.time << ","
							 << row.av_vm_stress << "," << row.av_vm_plastic_strain << "\n";

	   data_tree.put("Data.creep.macro_evolution", macro_evolution_csv.str());
	}

	// -------- aqs history ----------
	{
	   std::ostringstream plastic_events_csv;
	   plastic_events_csv << "index,element,eigenstrain_00,eigenstrain_11,eigenstrain_01\n";
	   for(auto &row : aqs_history.plastic)
		  plastic_events_csv << row.index << "," << row.element << "," << row.eigenstrain_00 << ","
							 << row.eigenstrain_11 << "," << row.eigenstrain_01 << "\n";

	   data_tree.put("Data.AQS.plastic_events", plastic_events_csv.str());


	   std::ostringstream driving_events_csv;
	   driving_events_csv << "index,dtime,dext_stress,dtotal_strain\n";
	   for(auto &row : aqs_history.driving)
		  driving_events_csv << row.index << "," << row.dtime << ","
							 << row.dext_stress << "," << row.dtotal_strain << "\n";

	   data_tree.put("Data.AQS.driving_events", driving_events_csv.str());


	   std::ostringstream macro_evolution_csv;
	   macro_evolution_csv << "index,ext_stress,total_strain,time,av_vm_stress,av_vm_plastic_strain\n";
	   for(auto &row : aqs_history.macro_evolution)
		  macro_evolution_csv << row.index << "," << row.ext_stress << ","
							 << row.total_strain << "," << row.time << ","
							 << row.av_vm_stress << "," << row.av_vm_plastic_strain << "\n";

	   data_tree.put("Data.AQS.macro_evolution", macro_evolution_csv.str());
	}


	// create a descriptive filename
   	std::ostringstream filename;

	filename << "Nx_" << p.Nx
	         << "+gamma_" << std::fixed << std::setprecision(2) << p.gamma
	         << "+lambda_" << p.lambda
	         << "+k_" << p.k
		     << "+G_" << p.G
		     << "+T_" << p.temperature
		     << "+ext_stress_" << p.ext_stress_creep
			 << "+seed_" << p.seed
			 << ".json";

   std::ofstream output_file( filename.str() );
   boost::property_tree::json_parser::write_json(output_file, data_tree);
   output_file.close();
}


void run(const Parameters &p, dealii::ConditionalOStream & cout)
{
   //----- SET UP ------

   constexpr unsigned dim = 2;
   std::mt19937 generator(p.seed);

   dealii::SymmetricTensor<4, dim> C = mepls::utils::tensor::make_isotropic_stiffness<dim>(p.G, p.nu);

   mepls::element::Vector<dim> elements;

   for(double n = 0; n < p.Nx * p.Ny; ++n)
   {
      example::element::Scalar<dim>::Config conf;
      conf.number = n;
      conf.gamma = p.gamma;
      conf.lambda = p.lambda;
      conf.k = p.k;
      conf.temperature = p.temperature;

      auto element = new example::element::Scalar<dim>(conf, generator);
      element->C(C);

      elements.push_back(element);
   }

   mepls::elasticity_solver::LeesEdwards<dim> solver(p.Nx, p.Ny, mepls::elasticity_solver::ControlMode::stress);
   for(auto &element : elements)
      solver.set_elastic_properties(element->number(), element->C());
   solver.setup_and_assembly();

   mepls::element::calculate_local_stress_coefficients_central(elements, solver);
   mepls::element::calculate_ext_stress_coefficients(elements, solver);

   mepls::system::Standard<dim> system(elements, solver, generator);


   //----- CREEP DEFORMATION ------

   mepls::history::History<dim> creep_history("creep_history");
   system.set_history(creep_history);
   creep_history.add_macro(system);

   // apply an external (stress) load of amplitude p.ext_stress_creep
   mepls::dynamics::fixed_load_increment(p.ext_stress_creep, system);
   creep_history.add_macro(system);

   mepls::dynamics::KMC<dim> kmc;

   mepls::utils::ContinueSimulation continue_creep;

	while( continue_creep() )
	{
		cout << system.macrostate["time"] << " " << system.macrostate["total_strain"] <<" "
			  << system.macrostate["ext_stress"] << std::endl;

		kmc(system);
		creep_history.add_macro(system);

		mepls::dynamics::relaxation(system, continue_creep);
		creep_history.add_macro(system);

	  continue_creep(system.macrostate["time"] < p.time_limit_creep, "creep time limit reached");
	}

	 cout << continue_creep << std::endl;


   //----- TRANSITION TO AQS ------

	// remove the external load
    mepls::dynamics::fixed_load_increment(-p.ext_stress_creep, system);
	creep_history.add_macro(system);

	// relax possible unstable slip system after the load change
	mepls::dynamics::relaxation(system, continue_creep);
	creep_history.add_macro(system);

    // in the AQS, we start measuring the macroscale strain from zero
	system.macrostate.clear();

	for(auto & element : elements)
		element->state_to_prestress();

	// we switch the driving mode to strain-controlled during AQS
    solver.set_control_mode(mepls::elasticity_solver::ControlMode::strain);

    // since the type of driving conditions have changed, the local stress change induced by a unit
    // load increment is different. We need to re-compute the ext_stress_coefficients
    mepls::element::calculate_ext_stress_coefficients(elements, solver);


    //----- ATHERMAL QUASISTATIC SHEAR ------

   	mepls::history::History<dim> aqs_history("aqs_history");
  	system.set_history(aqs_history);
   	aqs_history.add_macro(system);

   mepls::utils::ContinueSimulation continue_AQS_simulation;
   while(continue_AQS_simulation())
   {
		cout << system.macrostate["total_strain"] << " " << system.macrostate["ext_stress"] << std::endl;

      mepls::dynamics::extremal_dynamics_step(system);
      aqs_history.add_macro(system);

      mepls::dynamics::relaxation(system, continue_AQS_simulation);
      aqs_history.add_macro(system);

      continue_AQS_simulation(system.macrostate["total_strain"] < p.strain_limit_aqs, "AQS strain limit reached");
   }

	   cout << continue_AQS_simulation << std::endl;

   for(auto &element : elements)
      delete element;

	// write the simulation data, using its own dedicated function
	write_data(creep_history, aqs_history, p);
}


int main(int argc, char *argv[])
{
   // Read the command line arguments. We define the -f 'filename' to pass the
   // path to the parameters file
   cli::Parser parser(argc, argv);
   parser.set_optional<std::string>("f", "file", "./default.prm", "Name of the input configuration file");
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
      std::cout << "Configuration file " << parser.get<std::string>("f") << " created" << std::endl;
      return 1;
   }



	unsigned int n_rep = p.n_rep;
	if(n_rep < omp_get_max_threads())
		n_rep = omp_get_max_threads();

	// initialize the master engine with the master seed
	std::srand(p.seed);

	#pragma omp parallel
	{
		unsigned int n_threads = omp_get_max_threads();
		unsigned int id = omp_get_thread_num();
		unsigned int rep_per_thread = int( n_rep / n_threads );

		dealii::ConditionalOStream cout(std::cout, id==0 and p.verbose);

		for(unsigned int n = 0; n < rep_per_thread; ++n)
		{
			Parameters p_thread = p;
			#pragma critical
			{
				// generate a seed for a simulation run
				p_thread.seed = std::rand();
			};

			run(p_thread, cout);
		}

	}

   return 0;
}
