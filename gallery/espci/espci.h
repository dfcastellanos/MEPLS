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

#ifndef _ESPCI_H
#define _ESPCI_H

#include <mepls/patches.h>
#include <mepls/history.h>
#include <mepls/utils.h>
#include <mepls/solver.h>
#include <mepls/dynamics.h>
#include <mepls/system.h>
#include <mepls/element.h>
#include <mepls/snapshot.h>

#include <H5Cpp.h>

#include <omp.h>
#include <fstream>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/grid/grid_out.h>

/*! This namespace contains classes and functions to perform simulations based on ESPCI models. */
namespace espci
{

namespace parameters
{

struct Material
{
	// preeigenstrain
	bool init_eigenstrain = true;
	double init_eigenstrain_std_dev = 0.117;
	double init_eigenstrain_std_vol = 0.1;
	double init_eigenstrain_av_vol = 0.1;
	double init_eigenstrain_std_av_vol = 0.2;

	// threshold distribution
	double lambda = 5.;
	double lambda_quench = 5.;
	double k = 1.9;
	double k_quench = 1.9;
	double alpha_tau = 0.;

	// modulus distribution
	double G = 1.;
	double G_quench = 1.;
	double K = 1.;
	double K_quench = 1.;
	double gamma_pl_trans = 0.1;
	double beta = 0;

	// others
	double coupling_constant = 2.;
	double temperature_liquid = 0.2;
	double temperature_relaxation = 0.2;
	double activation_rate = 1.;
	unsigned int n_slip_systems = 1;

	double av_U = 0;
	double std_U = 1;
	double A = 1;
	double B = 1;

	void declare_entries(dealii::ParameterHandler &prm)
	{
		/*! Declare the variables of this struct into the input parser object. */

		prm.enter_subsection("Material");

		prm.declare_entry("init_eigenstrain_std_dev", mepls::utils::str::to_string(init_eigenstrain_std_dev),
						  dealii::Patterns::Double(), "");
		prm.declare_entry("init_eigenstrain_std_vol",
						  mepls::utils::str::to_string(init_eigenstrain_std_vol),
						  dealii::Patterns::Double(), "");
		prm.declare_entry("init_eigenstrain_av_vol",
						  mepls::utils::str::to_string(init_eigenstrain_av_vol),
						  dealii::Patterns::Double(), "");
		prm.declare_entry("init_eigenstrain_std_av_vol",
						  mepls::utils::str::to_string(init_eigenstrain_std_av_vol),
						  dealii::Patterns::Double(), "");

		prm.declare_entry("G", mepls::utils::str::to_string(G),
						  dealii::Patterns::Double(), "");
		prm.declare_entry("G_quench", mepls::utils::str::to_string(G_quench),
						  dealii::Patterns::Double(), "");
		prm.declare_entry("K", mepls::utils::str::to_string(K),
						  dealii::Patterns::Double(), "");
		prm.declare_entry("K_quench", mepls::utils::str::to_string(K_quench),
						  dealii::Patterns::Double(), "");
		prm.declare_entry("gamma_pl_trans", mepls::utils::str::to_string(gamma_pl_trans),
						  dealii::Patterns::Double(), "");
		prm.declare_entry("beta", mepls::utils::str::to_string(beta),
						  dealii::Patterns::Double(), "");

		prm.declare_entry("lambda", mepls::utils::str::to_string(lambda),
						  dealii::Patterns::Double(), "");
		prm.declare_entry("lambda_quench", mepls::utils::str::to_string(lambda_quench),
						  dealii::Patterns::Double(), "");
		prm.declare_entry("k", mepls::utils::str::to_string(k), dealii::Patterns::Double(), "");
		prm.declare_entry("k_quench", mepls::utils::str::to_string(k_quench),
						  dealii::Patterns::Double(), "");
		prm.declare_entry("alpha_tau", mepls::utils::str::to_string(alpha_tau),
						  dealii::Patterns::Double(0.0), "");
		prm.declare_entry("coupling_constant", mepls::utils::str::to_string(coupling_constant),
						  dealii::Patterns::Double(), "");

		prm.declare_entry("temperature_liquid", mepls::utils::str::to_string(temperature_liquid),
						  dealii::Patterns::Double(0.0), "");
		prm.declare_entry("temperature_relaxation", mepls::utils::str::to_string
		(temperature_relaxation),
						  dealii::Patterns::Double(0.0), "");
		prm.declare_entry("activation_rate", mepls::utils::str::to_string(activation_rate),
						  dealii::Patterns::Double(0.0), "");

		prm.declare_entry("n_slip_systems", mepls::utils::str::to_string(n_slip_systems),
						  dealii::Patterns::Integer(0), "");
		prm.declare_entry("init_eigenstrain", mepls::utils::str::to_string(init_eigenstrain),
						  dealii::Patterns::Bool(), "");

		prm.declare_entry("av_U", mepls::utils::str::to_string(av_U),
						  dealii::Patterns::Double(), "");
		prm.declare_entry("std_U", mepls::utils::str::to_string(std_U),
						  dealii::Patterns::Double(), "");
		prm.declare_entry("A", mepls::utils::str::to_string(A),
						  dealii::Patterns::Double(), "");
		prm.declare_entry("B", mepls::utils::str::to_string(B),
						  dealii::Patterns::Double(), "");

		prm.leave_subsection();

	}

	void load(dealii::ParameterHandler &prm)
	{
		/*! Load values parsed by the input parser object into the structs variables. */

		prm.enter_subsection("Material");

		init_eigenstrain_std_dev = prm.get_double("init_eigenstrain_std_dev");
		init_eigenstrain_std_vol = prm.get_double("init_eigenstrain_std_vol");
		init_eigenstrain_av_vol = prm.get_double("init_eigenstrain_av_vol");
		init_eigenstrain_std_av_vol = prm.get_double("init_eigenstrain_std_av_vol");

		G = prm.get_double("G");
		G_quench = prm.get_double("G_quench");
		K = prm.get_double("K");
		K_quench = prm.get_double("K_quench");
		gamma_pl_trans = prm.get_double("gamma_pl_trans");
		beta = prm.get_double("beta");

		lambda = prm.get_double("lambda");
		lambda_quench = prm.get_double("lambda_quench");
		k = prm.get_double("k");
		k_quench = prm.get_double("k_quench");
		alpha_tau = prm.get_double("alpha_tau");
		coupling_constant = prm.get_double("coupling_constant");
		temperature_liquid = prm.get_double("temperature_liquid");
		temperature_relaxation = prm.get_double("temperature_relaxation");
		activation_rate = prm.get_double("activation_rate");

		n_slip_systems = prm.get_integer("n_slip_systems");
		init_eigenstrain = prm.get_bool("init_eigenstrain");

		av_U = prm.get_double("av_U");
		std_U = prm.get_double("std_U");
		A = prm.get_double("A");
		B = prm.get_double("B");

		prm.leave_subsection();
	}

};



/*! Struct with parameters related to the simulation setup. */
struct Simulation
{
	unsigned int n_theta = 1;
	std::vector<int> N_patch_list;
	unsigned int n_rep = 1;

	unsigned int Nx = 32;
	/*!<  Number of elements in the x-direction of the rectangular lattice forming the system. */

	unsigned int Ny = 32;
	/*!< Number of elements in the y-direction of the rectangular lattice forming the system.  */

	unsigned int seed = 1234567;
	/*!<  Seed to initialized the random number engine. */

	std::string monitor_name = "total_strain";
	/*!< Magnitude used to check whether the simulation should stop and whether a snapshot should be
	 *  taken (see @ref snapshot). The posibilities correspond to the keys in
	 *  @ref MacroState.monitor_map. */

	double monitor_limit = 3.0;
	/*!<  Value of the magnitude defined by @ref monitor_name at which the simulation must stop. */

	bool parent_liquid = true;
	bool thermal_relaxation = true;
	bool reload = true;
	bool het_elasticity = false;
	bool do_ee = true;

	void declare_entries(dealii::ParameterHandler &prm)
	{
		/*! Declare the variables of this struct into the input parser object. */
		prm.enter_subsection("Simulation Setup");

		prm.declare_entry("n_rep", mepls::utils::str::to_string(n_rep),
						  dealii::Patterns::Integer(0), "");
		prm.declare_entry("n_theta", mepls::utils::str::to_string(n_theta),
						  dealii::Patterns::Integer(0), "");
		prm.declare_entry("N_patch_list", "", dealii::Patterns::FileName(), "");
		prm.declare_entry("Nx", mepls::utils::str::to_string(Nx), dealii::Patterns::Integer(0), "");
		prm.declare_entry("Ny", mepls::utils::str::to_string(Ny), dealii::Patterns::Integer(0), "");
		prm.declare_entry("seed", mepls::utils::str::to_string(seed), dealii::Patterns::Integer(0),
						  "");
		prm.declare_entry("monitor_name", monitor_name, dealii::Patterns::Selection(
			"av_vm_plastic_strain|load|ext_stress|time|total_strain"), "");
		prm.declare_entry("monitor_limit", mepls::utils::str::to_string(monitor_limit),
						  dealii::Patterns::Double(0.0), "");
		prm.declare_entry("parent_liquid", mepls::utils::str::to_string(parent_liquid),
						  dealii::Patterns::Bool(), "");
		prm.declare_entry("thermal_relaxation", mepls::utils::str::to_string(thermal_relaxation),
						  dealii::Patterns::Bool(), "");
		prm.declare_entry("reload", mepls::utils::str::to_string(reload), dealii::Patterns::Bool(),
						  "");
		prm.declare_entry("het_elasticity", mepls::utils::str::to_string(het_elasticity),
						  dealii::Patterns::Bool(), "");
		prm.declare_entry("do_ee", mepls::utils::str::to_string(do_ee), dealii::Patterns::Bool(),
						  "");

		prm.leave_subsection();
	}

	void load(dealii::ParameterHandler &prm)
	{
		/*! Load values parsed by the input parser object into the structs variables. */
		prm.enter_subsection("Simulation Setup");

		n_rep = prm.get_integer("n_rep");
		n_theta = prm.get_integer("n_theta");
		Nx = prm.get_integer("Nx");
		Ny = prm.get_integer("Ny");
		seed = prm.get_integer("seed");
		monitor_name = prm.get("monitor_name");
		monitor_limit = prm.get_double("monitor_limit");
		N_patch_list = mepls::utils::str::parse_list_integers(prm.get("N_patch_list"));
		parent_liquid = prm.get_bool("parent_liquid");
		thermal_relaxation = prm.get_bool("thermal_relaxation");
		reload = prm.get_bool("reload");
		het_elasticity = prm.get_bool("het_elasticity");
		do_ee = prm.get_bool("do_ee");

		prm.leave_subsection();
	}
};


/*! Struct with parameters related to the simulation output. */
struct Output
{
	std::string filename = "";
	/*!< Name of the output file. If the name is auto, then the name will be automatically generated
	 * based on the simulation parameters. */

	bool verbosity = true;
	/*!< If true, show output on the screen during an on-going simulation. */

	std::string path = "./results";
	/*!< Path to the directory where the output files should be written. */

	std::string snapshots = "";
	/*!< Select the type of snapshots to take. The string must contains some of the following types:
	 *  threshold, stress, def_grad, patches. */

	double snapshots_min = 0.;
	/*!< Minimum value of the magnitude @ref Simulation.monitor_name from which snapshots should
	 * start to be taken (see @ref snapshot).  */

	double snapshots_max = 1.;
	/*!<  Maximum value of the magnitude @ref Simulation.monitor_name until which snapshots should
	 * be taken (see @ref snapshot). */

	double snapshots_interval = 0.05;
	/*!< Interval to create a list of value from @ref snapshots_min to @ref snapshots_max at which
	 *  snapshots should be taken (see @ref snapshot).   */

	double snapshots_sensitivity = 0.01;

	/*!< Numerical tolerance to decide whether the value of @ref Simulation.monitor_name is that at
	 * which a snapshot must be taken (see @ref snapshot). */

	void declare_entries(dealii::ParameterHandler &prm)
	{
		/*! Declare the variables of this struct into the input parser object. */
		prm.enter_subsection("Output");

		prm.declare_entry("filename", filename, dealii::Patterns::FileName(), "");
		prm.declare_entry("verbosity", mepls::utils::str::to_string(verbosity),
						  dealii::Patterns::Bool(), "");
		prm.declare_entry("output_path", path, dealii::Patterns::FileName(), "");
		prm.declare_entry("snapshots", snapshots, dealii::Patterns::FileName(), "");
		prm.declare_entry("snapshots_min", mepls::utils::str::to_string(snapshots_min),
						  dealii::Patterns::Double(0.0), "");
		prm.declare_entry("snapshots_max", mepls::utils::str::to_string(snapshots_max),
						  dealii::Patterns::Double(0.0), "");
		prm.declare_entry("snapshots_interval", mepls::utils::str::to_string(snapshots_interval),
						  dealii::Patterns::Double(0.0), "");
		prm.declare_entry("snapshots_sensitivity",
						  mepls::utils::str::to_string(snapshots_sensitivity),
						  dealii::Patterns::Double(0.0), "");

		prm.leave_subsection();
	}

	void load(dealii::ParameterHandler &prm)
	{
		/*! Load values parsed by the input parser object into the structs variables. */
		prm.enter_subsection("Output");

		filename = prm.get("filename");
		verbosity = prm.get_bool("verbosity");
		path = prm.get("output_path");
		snapshots = prm.get("snapshots");
		snapshots_min = prm.get_double("snapshots_min");
		snapshots_max = prm.get_double("snapshots_max");
		snapshots_interval = prm.get_double("snapshots_interval");
		snapshots_sensitivity = prm.get_double("snapshots_sensitivity");

		prm.leave_subsection();
	}
};


/*! A struct with all parameters necessary for running a simulation. It is composed of other structs
 *  defining the parameters of specific parts such as material parameters, simulation setup
 *  parameters, and output setup parameters. */
struct Parameters
{
	void load_file(const std::string &filename)
	{
		/*! Create a parser object and forward it to the structs so they load their own
		 * variables. */

		dealii::ParameterHandler prm;

		sim.declare_entries(prm);
		out.declare_entries(prm);
		mat.declare_entries(prm);

		prm.parse_input(filename);

		sim.load(prm);
		out.load(prm);
		mat.load(prm);
	}

	void generate_file(const std::string &filename)
	{
		/*! Write a text file with the right structure, such that it can be parsed for loading this
		 *  struct. */

		std::ofstream outfile(filename);
		dealii::ParameterHandler prm;

		sim.declare_entries(prm);
		out.declare_entries(prm);
		mat.declare_entries(prm);

		prm.print_parameters(outfile, dealii::ParameterHandler::OutputStyle::Text);
	}

	Simulation sim;
	/*!< Struct with parameters related to the simulation setup. */

	Output out;
	/*!< Struct with parameters related to the simulation output. */

	Material mat;
	/*!< Struct with parameters related to material properties. */

};

} // parameters


namespace slip
{
/*! This class represents a slip system defined by a slip angle and a slip threshold.
*
* \note In contrast with usual crystal slip systems, this class considers only positive slip, i.e.,
 * if the angle is 0 rad the slip plane is the 0.3 rad plane; the negative version of that slip
 * would correspond to the activation of another slip system (represented by a different object)
 * oriented with 0.3+pi/2 rad. */
template<int dim>
class Oriented: public mepls::slip::Slip<dim>
{
public:

	struct Config
	{
		double threshold = 0.;
		double angle = 0.;
		double coupling_constant = 2.;
		double alpha_tau = 0.;
		double activation_rate_0 = 1.;
		double temperature = 0.2;
	};

	Oriented(std::mt19937 &generator_, std::uniform_real_distribution<double> &unif_dist_, const
	Config &conf_)
		:
		mepls::slip::Slip<dim>(),
		conf(conf_),
		threshold_0(conf_.threshold),
		unif_distribution(unif_dist_),
		generator(generator_)
	{
		M_Assert(threshold_0 > 0., "Expected threshold > 0");

		angle = mepls::utils::mod(conf.angle, M_PI);
		M = mepls::utils::tensor::make_schmid<dim>(angle);M_Assert(M[0][1] == M[1][0], "");M_Assert(
			M[0][0] == -M[1][1], "");
	};

	Oriented(Oriented *slip)
		:
		mepls::slip::Slip<dim>(slip),
		M(slip->M),
		conf(slip->conf),
		threshold_0(slip->threshold_0),
		unif_distribution(slip->unif_distribution),
		generator(slip->generator)
	{
	}

	void update() override
	{
		eff_shear_stress = M * parent->stress();
		pressure = dealii::trace(parent->stress()) / double(dim);

		threshold = threshold_0 - conf.alpha_tau * pressure;
		if(threshold < 0.)
			threshold = 1e-3;

		barrier = threshold - eff_shear_stress;
		activation_rate = conf.activation_rate_0 * std::exp(-barrier / conf.temperature);

		M_Assert(threshold > 0., "Expected threshold > 0");
		M_Assert(not std::isnan(threshold),	"Threshold is NaN");
		M_Assert(not std::isnan(eff_shear_stress), "eff_shear_stress is NaN");
	}

	double get_critical_load_increment() override
	{
		M_Assert(barrier >= 0, "The barrier cannot be negative when computing the critical load increment.");

		return barrier / (M * parent->ext_stress_coeff() + 0.5 * conf.alpha_tau * dealii::trace(
			parent->ext_stress_coeff()));
	}

	dealii::SymmetricTensor<2, dim> get_eigenstrain_increment() override
	{
		bool thermal = conf.temperature > 0;
		double eff_shear_stress_variation = get_local_shear_stress_variation();
		double A = M * parent->S() * M;
		double gamma = 2 * eff_shear_stress_variation / A;
		dealii::SymmetricTensor<2, dim> eigenstrain_dev = gamma * M;

		return eigenstrain_dev;
	}

	double get_local_shear_stress_variation()
	{
		bool thermal = conf.temperature > 0;
		double eff_shear_stress_variation = 0.;
		
		std::gamma_distribution<double> g_a(1, 1);
		std::gamma_distribution<double> g_b(conf.coupling_constant, 1);
		double za = g_a(generator);
		double zb = g_b(generator);

		if(thermal)
		{
			eff_shear_stress_variation = -(eff_shear_stress + thermal * (barrier + conf.temperature)) * za / (za + zb);
		}else{
			eff_shear_stress_variation = - eff_shear_stress * za/(za+zb); // c=-1.254 for stat fit
			// eff_shear_stress_variation = - threshold * za/(za+zb); // c=-0.812 for stat fit

			M_Assert(eff_shear_stress *eff_shear_stress_variation<1,
					 "Eff. shear stress variation does not reduce local stress. Is the orientation of plastic event according to local stress?");
		}
		return eff_shear_stress_variation;
	}

	Oriented<dim> *make_copy_impl() override
	{
		return new Oriented<dim>(this);
	}

	Oriented<dim> *make_copy()
	{
		return make_copy_impl();
	}

	using mepls::slip::Slip<dim>::angle;
	using mepls::slip::Slip<dim>::eff_shear_stress;
	using mepls::slip::Slip<dim>::pressure;
	using mepls::slip::Slip<dim>::threshold;
	using mepls::slip::Slip<dim>::barrier;
	using mepls::slip::Slip<dim>::parent;
	using mepls::slip::Slip<dim>::activation_rate;

	dealii::SymmetricTensor<2, dim> M;
	Config conf;
	double threshold_0;
	std::uniform_real_distribution<double> &unif_distribution;
	std::mt19937 &generator;
};

} // namespace slip


namespace element
{

template<int dim>
class Anisotropic: public mepls::element::Element<dim>
{
public:

	struct Config
	{
		double G = 1.;
		double K = 1.;
		double G_quench = 1.;
		double K_quench = 1.;

		double alpha_tau = 0.;
		double coupling_constant = 2.;
		double temperature = 0.;
		double activation_rate_0 = 1.;
		double k = 2;
		double k_quench = 2;
		double lambda = 1.;
		double lambda_quench = 1.;
		double gamma_pl_trans = 1e-5;
		double beta = 0;
		bool vary_elastic_properties = false;

		unsigned int n_slip_systems = 1;
		unsigned int number = 0;

		double av_U = 0;
		double std_U = 0;
		double A = 0;
		double B = 0;
	};


	Anisotropic(std::mt19937 &generator_,
				const Config &conf_)
		:
		mepls::element::Element<dim>(),
		generator(generator_),
		unif_distribution(0, 1),
		normal_distribution(0,1),
		conf(conf_)
	{
		this->number(conf.number);

		mepls::event::Plastic<dim> plastic_event;
		renew_thresholds(plastic_event);
	}

	void renew_structural_properties_impl(mepls::event::Plastic<dim> &plastic_event) override
	{
		if(plastic_event.renew_slip_properties)
			renew_thresholds(plastic_event);
	}

	void renew_thresholds(mepls::event::Plastic<dim> &plastic_event)
	{
		this->remove_slip_systems();

		double k, lambda, G, K;

		if(this->integrated_vm_eigenstrain()==0. or not conf.vary_elastic_properties)
		{
			lambda = conf.lambda_quench;
			k = conf.k_quench;
			G = conf.G_quench;
			K = conf.K_quench;
		}
		else
		{
			assert(plastic_event.dplastic_strain > 0.);
			double x = std::exp(-(plastic_event.dplastic_strain/(conf.gamma_pl_trans*std::pow(lambda_old/conf.lambda,conf.beta))));
			lambda = (lambda_old - conf.lambda)*x + conf.lambda;
			k = (k_old - conf.k)*x + conf.k;
			G = (G_old - conf.G)*x + conf.G;
			K = (K_old - conf.K)*x + conf.K;
		}

		auto C = mepls::utils::tensor::make_mandel_anisotropic_stiffness<dim>(K, G, G, 0.);
		this->C(mepls::utils::tensor::mandel_to_standard_rank4<dim>(C));

		typename slip::Oriented<dim>::Config slip_conf;
		slip_conf.alpha_tau = conf.alpha_tau;
		slip_conf.coupling_constant = conf.coupling_constant;
		slip_conf.activation_rate_0 = conf.activation_rate_0;
		slip_conf.temperature = conf.temperature;

//		slip_conf.angle = 0;
//		slip_conf.threshold = mepls::utils::rand::get_weibull_rand(k, lambda, unif_distribution(generator));
//		this->add_slip_system(new slip::Oriented<dim>(generator, unif_distribution, slip_conf));
//
//		slip_conf.angle = M_PI / 2.;
//		this->add_slip_system(new slip::Oriented<dim>(generator, unif_distribution, slip_conf));

		for(unsigned int i = 0; i < conf.n_slip_systems; ++i)
		{
			double alpha2 = unif_distribution(generator) * M_PI;

			slip_conf.angle = alpha2;
			slip_conf.threshold = mepls::utils::rand::get_weibull_rand(k, lambda, unif_distribution(generator));
			this->add_slip_system(new slip::Oriented<dim>(generator, unif_distribution, slip_conf));

			slip_conf.angle = alpha2 + M_PI / 2.;
			slip_conf.threshold = mepls::utils::rand::get_weibull_rand(k, lambda, unif_distribution(generator));
			this->add_slip_system(new slip::Oriented<dim>(generator, unif_distribution, slip_conf));

			slip_conf.angle = alpha2 + M_PI / 4.;
			slip_conf.threshold = mepls::utils::rand::get_weibull_rand(k, lambda, unif_distribution(generator));
			this->add_slip_system(new slip::Oriented<dim>(generator, unif_distribution, slip_conf));

			slip_conf.angle = alpha2 + M_PI / 4. + M_PI / 2.;
			slip_conf.threshold = mepls::utils::rand::get_weibull_rand(k, lambda, unif_distribution(generator));
			this->add_slip_system(new slip::Oriented<dim>(generator, unif_distribution, slip_conf));
		}

		lambda_old = lambda;
		k_old = k;
		G_old = G;
		K_old = K;
	}

	Anisotropic<dim> *make_copy_impl() override
	{
		auto new_element = new Anisotropic<dim>(this->generator, this->conf);

		new_element->make_copy(this);

		return new_element;
	}

	Anisotropic<dim> *make_copy()
	{
		return make_copy_impl();
	}

	void make_copy(Anisotropic<dim> *input_element)
	{
		conf = input_element->conf;
		mepls::element::Element<dim>::make_copy(input_element);
	}

	const Config &config() const
	{
		return conf;
	}

	void config(const Config &conf_)
	{
		conf = conf_;

		// slip systems must be informed about changes in the element configuration (slips have
		// pointer to parent element::Element so they cannot access element.conf
		// directly) this
		// operation is safe as long as we know that Anisotropic element type has Oriented slips
		// (which is always the case)
		for(auto &slip : *this)
		{
			auto &slip_conf = static_cast<slip::Oriented<dim> *>(slip)->conf;
			slip_conf.alpha_tau = conf.alpha_tau;
			slip_conf.coupling_constant = conf.coupling_constant;
			slip_conf.temperature = conf.temperature;
			slip_conf.activation_rate_0 = conf.activation_rate_0;
		}

		// slips must be informed. For example, if the temperature has changed,
		// each slip must update its activation rate
		for(auto &slip : *this)
			slip->update();
	}

protected:
	std::mt19937 &generator;
	std::uniform_real_distribution<double> unif_distribution;
	std::normal_distribution<double> normal_distribution;
	Config conf;
	double lambda_old, k_old, G_old, K_old;
};

} // namespace element


namespace history
{

struct MacroSummaryRow : public mepls::history::MacroSummaryRow
{
	/*! Struct to store system-scale properties */
	double G = 0.;
	double K = 0.;
};

template<int dim>
class History : public mepls::history::History<dim>
{
  public:

	History(parameters::Parameters p_, const std::string &inputname_ = "history")
		:
	mepls::history::History<dim>(inputname_),
	p(p_)
	{
		/*! Constructor. */
	}


	void add_macro(const mepls::system::System<dim> &system) override
	{
		if(mepls::history::History<dim>::closed_)
			return;

		MacroSummaryRow data;

		dealii::SymmetricTensor<4,dim> av_C;
		for(auto &element : system)
			av_C += element->C();
		av_C /= double(system.size());

		data.G = av_C[0][1][0][1];
		data.K = (av_C[0][0][0][0]+av_C[0][0][1][1])/2.;


		mepls::history::History<dim>::add_macro_default(system, data);

		macro_evolution_espci.push_back(data);
	}

	std::vector<MacroSummaryRow> macro_evolution_espci;

	parameters::Parameters p;
};


} // history

template<int dim>
class GlobalPropertiesSnapshot
{
public:

	struct DataRow
	{
		int n = 0; // sample number

		float C_00 = 0.; // stiffness tensor
		float C_01 = 0.;
		float C_02 = 0.;
		float C_10 = 0.;
		float C_11 = 0.;
		float C_12 = 0.;
		float C_20 = 0.;
		float C_21 = 0.;
		float C_22 = 0.;

		float theta = 0.; // first instability
		float oi_eps = 0.;
		float ss_00 = 0.;
		float ss_11 = 0.;
		float ss_01 = 0.;
		float oi_00 = 0.;
		float oi_11 = 0.;
		float oi_01 = 0.;
	};

	GlobalPropertiesSnapshot(mepls::system::System<dim> &system,
							 mepls::elasticity_solver::LeesEdwards<dim> &solver,
							 unsigned int output_index_,
							 std::string monitor_mag_,
							 double desired_target_,
							 double recorded_target_)
		:
		recorded_mag("global_properties"),
		monitor_name(monitor_mag_),
		desired_target(desired_target_),
		recorded_target(recorded_target_),
		output_index(output_index_)
	{
		DataRow global_properties;

		// save full-system stiffness tensor
		double ext_gamma = 1.;
		solver.add_load_increment(ext_gamma / 2.);
		solver.solve();M_Assert(ext_gamma == 2. * solver.get_total_strain(), "");

		auto C = mepls::utils::tensor::standard_rank4_to_voigt<dim>(solver.get_global_stiffness());
		global_properties.C_00 = C[0][0];
		global_properties.C_01 = C[0][1];
		global_properties.C_02 = C[0][2];
		global_properties.C_10 = C[1][0];
		global_properties.C_11 = C[1][1];
		global_properties.C_12 = C[1][2];
		global_properties.C_20 = C[2][0];
		global_properties.C_21 = C[2][1];
		global_properties.C_22 = C[2][2];

		solver.clear();


		// compute full-system first instability the following will change the fields in the
		// elements. Therefore, we work on copies so the originals can be used later by other
		// functions
		mepls::element::Vector<dim> copy_elements;
		for(auto &element : system)
			copy_elements.push_back(element->make_copy());

		mepls::utils::ContinueSimulation continue_shear_test;

		mepls::element::calculate_ext_stress_coefficients(copy_elements, solver);

		auto system_replica = system.get_new_instance(copy_elements, solver, system.generator);

		// we can use the tools for patches to study the full system response
		mepls::patches::PatchPropertiesTensorial<dim> patch_properties;

		// external loading is shear with only epsxy!=0.
		double theta = 0.;

		// to compare with the MD full system global_properties, the external strain discrete
		// incremement is here 1e-4 (different from the 1e-3 used for the patches)
		double dgamma = 1e-4;
		mepls::patches::apply_patch_shear_test<dim>(patch_properties, *system_replica,
													continue_shear_test, false, dgamma);

		global_properties.oi_eps = patch_properties.resolved_elastic_shear_strain_oi;
		global_properties.ss_00 = patch_properties.stress_ss[0][0];
		global_properties.ss_11 = patch_properties.stress_ss[1][1];
		global_properties.ss_01 = patch_properties.stress_ss[0][1];
		global_properties.oi_00 = patch_properties.stress_oi[0][0];
		global_properties.oi_11 = patch_properties.stress_oi[1][1];
		global_properties.oi_01 = patch_properties.stress_oi[0][1];

		solver.clear();


		data.push_back(global_properties);
	};

	std::vector<DataRow> data;
	/*!< Container to store the recorde data. */

	std::string recorded_mag;
	/*!< Name of the field storaged in @ref data. */

	std::string monitor_name;
	/*!< Name of the magnitude used to check whether the snapshot should be taken or not. */

	double desired_target;
	/*!< Value of the @ref monitor_name at which we desired to take the snapshot. */

	double recorded_target;
	/*!< Value of the @ref monitor_name at which the snapshot is actually taken. */

	unsigned int output_index;
	/*!< Global event index from the @ref mepls::event::History at which the snapshot is taken. */
};


namespace write
{

inline void file_attrs(H5::H5File &file, const parameters::Parameters &p)
{
	H5::DataSpace att_space(H5S_SCALAR);
	file.createAttribute("lambda", H5::PredType::NATIVE_DOUBLE, att_space).write(
		H5::PredType::NATIVE_DOUBLE, &p.mat.lambda);
	file.createAttribute("k", H5::PredType::NATIVE_DOUBLE, att_space).write(
		H5::PredType::NATIVE_DOUBLE, &p.mat.k);
	file.createAttribute("lambda_quench", H5::PredType::NATIVE_DOUBLE, att_space).write(
		H5::PredType::NATIVE_DOUBLE, &p.mat.lambda_quench);
	file.createAttribute("k_quench", H5::PredType::NATIVE_DOUBLE, att_space).write(
		H5::PredType::NATIVE_DOUBLE, &p.mat.k_quench);
	file.createAttribute("gamma_pl_trans", H5::PredType::NATIVE_DOUBLE, att_space).write(
		H5::PredType::NATIVE_DOUBLE, &p.mat.gamma_pl_trans);
	file.createAttribute("beta", H5::PredType::NATIVE_DOUBLE, att_space).write(
		H5::PredType::NATIVE_DOUBLE, &p.mat.beta);
	file.createAttribute("alpha_tau", H5::PredType::NATIVE_DOUBLE, att_space).write(
		H5::PredType::NATIVE_DOUBLE, &p.mat.alpha_tau);
	file.createAttribute("G", H5::PredType::NATIVE_DOUBLE, att_space).write(
		H5::PredType::NATIVE_DOUBLE, &p.mat.G);
	file.createAttribute("K", H5::PredType::NATIVE_DOUBLE, att_space).write(
		H5::PredType::NATIVE_DOUBLE, &p.mat.K);
	file.createAttribute("G_quench", H5::PredType::NATIVE_DOUBLE, att_space).write(
		H5::PredType::NATIVE_DOUBLE, &p.mat.G_quench);
	file.createAttribute("K_quench", H5::PredType::NATIVE_DOUBLE, att_space).write(
		H5::PredType::NATIVE_DOUBLE, &p.mat.K_quench);
	file.createAttribute("Nx", H5::PredType::NATIVE_UINT, att_space).write(
		H5::PredType::NATIVE_UINT, &p.sim.Nx);
	file.createAttribute("Ny", H5::PredType::NATIVE_UINT, att_space).write(
		H5::PredType::NATIVE_UINT, &p.sim.Ny);
	file.createAttribute("seed", H5::PredType::NATIVE_UINT, att_space).write(
		H5::PredType::NATIVE_UINT, &p.sim.seed);
	file.createAttribute("coupling_constant", H5::PredType::NATIVE_DOUBLE, att_space).write(
		H5::PredType::NATIVE_DOUBLE, &p.mat.coupling_constant);
	file.createAttribute("temperature_liquid", H5::PredType::NATIVE_DOUBLE, att_space).write(
		H5::PredType::NATIVE_DOUBLE, &p.mat.temperature_liquid);
	file.createAttribute("temperature_relaxation", H5::PredType::NATIVE_DOUBLE, att_space).write(
		H5::PredType::NATIVE_DOUBLE, &p.mat.temperature_relaxation);
	file.createAttribute("activation_rate", H5::PredType::NATIVE_DOUBLE, att_space).write(
		H5::PredType::NATIVE_DOUBLE, &p.mat.activation_rate);
	file.createAttribute("init_eigenstrain_av_vol", H5::PredType::NATIVE_DOUBLE, att_space).write(
		H5::PredType::NATIVE_DOUBLE, &p.mat.init_eigenstrain_av_vol);
	file.createAttribute("init_eigenstrain_std_av_vol", H5::PredType::NATIVE_DOUBLE, att_space).write(
		H5::PredType::NATIVE_DOUBLE, &p.mat.init_eigenstrain_std_av_vol);
	file.createAttribute("init_eigenstrain_std_vol", H5::PredType::NATIVE_DOUBLE, att_space).write(
		H5::PredType::NATIVE_DOUBLE, &p.mat.init_eigenstrain_std_vol);
	file.createAttribute("init_eigenstrain_std_dev", H5::PredType::NATIVE_DOUBLE, att_space).write(
		H5::PredType::NATIVE_DOUBLE, &p.mat.init_eigenstrain_std_dev);
	file.createAttribute("n_slip_systems", H5::PredType::NATIVE_UINT, att_space).write(
		H5::PredType::NATIVE_UINT, &p.mat.n_slip_systems);

	{
		H5std_string strwritebuf;
		strwritebuf = p.sim.monitor_name;
		H5::StrType strdatatype(H5::PredType::C_S1, strwritebuf.size());
		file.createAttribute("monitor_name", strdatatype, att_space).write(strdatatype,
																		   strwritebuf);
	}
}


template<int dim>
inline void evolution_history(H5::H5File &file,
						  const typename history::History<dim> &history)
{
	std::string path = "/"+history.name();
	if(not H5Lexists(file.getId(), path.c_str(), H5P_DEFAULT))
		file.createGroup(path.c_str());

	{   /* --------- write macro evolution ----------- */

		using DataRow = history::MacroSummaryRow;

		H5::CompType mtype(sizeof(DataRow));
		mtype.insertMember("time", HOFFSET(DataRow, time), H5::PredType::NATIVE_DOUBLE);
		mtype.insertMember("total_strain", HOFFSET(DataRow, total_strain), H5::PredType::NATIVE_DOUBLE);
		mtype.insertMember("ext_stress", HOFFSET(DataRow, ext_stress), H5::PredType::NATIVE_DOUBLE);
		mtype.insertMember("av_vm_plastic_strain", HOFFSET(DataRow, av_vm_plastic_strain),
						   H5::PredType::NATIVE_DOUBLE);
		mtype.insertMember("std_vm_plastic_strain", HOFFSET(DataRow, std_vm_plastic_strain),
						   H5::PredType::NATIVE_DOUBLE);
		mtype.insertMember("av_vm_stress", HOFFSET(DataRow, av_vm_stress), H5::PredType::NATIVE_DOUBLE);
		mtype.insertMember("std_vm_stress", HOFFSET(DataRow, std_vm_stress),
						   H5::PredType::NATIVE_DOUBLE);
		mtype.insertMember("av_energy_el", HOFFSET(DataRow, av_energy_el),
						   H5::PredType::NATIVE_DOUBLE);
		mtype.insertMember("std_energy_el", HOFFSET(DataRow, std_energy_el),
						   H5::PredType::NATIVE_DOUBLE);
		mtype.insertMember("av_energy_conf", HOFFSET(DataRow, av_energy_conf),
						   H5::PredType::NATIVE_DOUBLE);
		mtype.insertMember("std_energy_conf", HOFFSET(DataRow, std_energy_conf),
						   H5::PredType::NATIVE_DOUBLE);
		mtype.insertMember("av_stress_00", HOFFSET(DataRow, av_stress_00),
						   H5::PredType::NATIVE_DOUBLE);
		mtype.insertMember("std_stress_00", HOFFSET(DataRow, std_stress_00),
						   H5::PredType::NATIVE_DOUBLE);
		mtype.insertMember("av_stress_11", HOFFSET(DataRow, av_stress_11),
						   H5::PredType::NATIVE_DOUBLE);
		mtype.insertMember("std_stress_11", HOFFSET(DataRow, std_stress_11),
						   H5::PredType::NATIVE_DOUBLE);
		mtype.insertMember("av_stress_01", HOFFSET(DataRow, av_stress_01),
						   H5::PredType::NATIVE_DOUBLE);
		mtype.insertMember("std_stress_01", HOFFSET(DataRow, std_stress_01),
						   H5::PredType::NATIVE_DOUBLE);
		mtype.insertMember("index", HOFFSET(DataRow, index), H5::PredType::NATIVE_UINT);
		mtype.insertMember("shear_modulus", HOFFSET(DataRow, G), H5::PredType::NATIVE_DOUBLE);
		mtype.insertMember("bulk_modulus", HOFFSET(DataRow, K), H5::PredType::NATIVE_DOUBLE);

		hsize_t d[] = {history.macro_evolution_espci.size()};
		H5::DataSpace space(1, d);
		H5::DataSet dataset = file.createDataSet(path + "/macro_evolution", mtype, space);

		dataset.write(history.macro_evolution_espci.data(), mtype);
	}

	{   /* --------- write driving event history ----------- */

		using DataRow = typename mepls::history::DrivingRow;
		H5::CompType mtype(sizeof(DataRow));
		mtype.insertMember("index", HOFFSET(DataRow, index), H5::PredType::NATIVE_UINT);
		mtype.insertMember("dload", HOFFSET(DataRow, dload), H5::PredType::NATIVE_DOUBLE);
		mtype.insertMember("dext_stress", HOFFSET(DataRow, dext_stress),
						   H5::PredType::NATIVE_DOUBLE);
		mtype.insertMember("dpressure", HOFFSET(DataRow, dpressure), H5::PredType::NATIVE_DOUBLE);
		mtype.insertMember("dtotal_strain", HOFFSET(DataRow, dtotal_strain),
						   H5::PredType::NATIVE_DOUBLE);
		mtype.insertMember("dtime", HOFFSET(DataRow, dtime), H5::PredType::NATIVE_DOUBLE);
		mtype.insertMember("activation_protocol", HOFFSET(DataRow, activation_protocol),
						   H5::PredType::NATIVE_UINT);

		hsize_t d[] = {history.driving.size()};
		H5::DataSpace space(1, d);
		H5::DataSet dataset = file.createDataSet(path + "/driving_events", mtype, space);

		dataset.write(history.driving.data(), mtype);
	}

	{   /* --------- write plastic event history ----------- */

		using DataRow = typename mepls::history::PlasticRow;
		H5::CompType mtype(sizeof(DataRow));
		mtype.insertMember("index", HOFFSET(DataRow, index), H5::PredType::NATIVE_UINT);
		mtype.insertMember("element", HOFFSET(DataRow, element), H5::PredType::NATIVE_UINT);
		mtype.insertMember("eigenstrain_00", HOFFSET(DataRow, eigenstrain_00),
						   H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("eigenstrain_11", HOFFSET(DataRow, eigenstrain_11),
						   H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("eigenstrain_01", HOFFSET(DataRow, eigenstrain_01),
						   H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("acting_stress_00", HOFFSET(DataRow, acting_stress_00),
						   H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("acting_stress_11", HOFFSET(DataRow, acting_stress_11),
						   H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("acting_stress_01", HOFFSET(DataRow, acting_stress_01),
						   H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("slip_threshold", HOFFSET(DataRow, slip_threshold),
						   H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("dplastic_strain", HOFFSET(DataRow, dplastic_strain),
						   H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("activation_protocol", HOFFSET(DataRow, activation_protocol),
						   H5::PredType::NATIVE_UINT);

		hsize_t d[] = {history.plastic.size()};
		H5::DataSpace space(1, d);
		H5::DataSet dataset = file.createDataSet(path + "/plastic_events", mtype, space);

		dataset.write(history.plastic.data(), mtype);
	}

	//	   {   /* --------- write renewal event history ----------- */
	//		  using DataRow = typename mepls::history::History<dim>::RenewSlipRow;
	//	      H5::CompType mtype( sizeof(DataRow) );
	//	      mtype.insertMember( "index", HOFFSET(DataRow, index), H5::PredType::NATIVE_UINT);
	//	      mtype.insertMember( "element", HOFFSET(DataRow, element), H5::PredType::NATIVE_UINT);
	//	      mtype.insertMember( "threshold", HOFFSET(DataRow, threshold), H5::PredType::NATIVE_FLOAT);
	//	      mtype.insertMember( "slip_angle", HOFFSET(DataRow, slip_angle), H5::PredType::NATIVE_FLOAT);
	//
	//	      hsize_t d[] = {history.renew.size()};
	//	      H5::DataSpace space( 1, d );
	//	      H5::DataSet dataset = file.createDataSet(path+"/renew", mtype, space);
	//
	//	      dataset.write( history.renew.data(), mtype );
	//	   }

}


template<int dim>
inline void patch_info(H5::H5File &file,
					   std::string path,
					   const std::map<unsigned int,
									  std::vector<std::vector<unsigned int>>> &patch_to_element_map)
{
	if(not H5Lexists(file.getId(), path.c_str(), H5P_DEFAULT))
		file.createGroup(path);

	// the dataset index corresponds to the patch ref. element number;
	// each column number corresponds to the element number within the patch numbering scheme;
	// each value corresponds to the number of the element in the whole-system numbering scheme;

	for(auto &x : patch_to_element_map)
	{
		unsigned int n_patch = x.first;
		std::vector<std::vector<unsigned int>> v = x.second;

		unsigned int nrow = v.size();
		unsigned int ncol = v[0].size();
		int varray[nrow][ncol];
		for(unsigned int i = 0; i < nrow; ++i)
			for(unsigned int j = 0; j < ncol; ++j)
				varray[i][j] = v[i][j];

		hsize_t dimsf[2];
		dimsf[0] = v.size();
		dimsf[1] = v[0].size();
		H5::DataSpace dataspace(2, dimsf);

		H5::DataType datatype(H5::PredType::NATIVE_UINT);
		H5::DataSet dataset = file.createDataSet("/patch_info/" + std::to_string(n_patch), datatype,
												 dataspace);

		dataset.write(varray, H5::PredType::NATIVE_UINT);
		dataset.close();
		dataspace.close();
	}

}


template<int dim>
inline void element_info(H5::H5File &file,
						 std::string path,
						 const mepls::element::Vector<dim> &elements)
{   /* --------- write element setup ----------- */
	using ElementSetupRow = mepls::element::SetupRow<dim>;

	std::vector<ElementSetupRow> element_setup;
	for(auto &element : elements)
	{
		ElementSetupRow setup;
		setup.prestress_00 = element->prestress()[0][0];
		setup.prestress_11 = element->prestress()[1][1];
		setup.prestress_01 = element->prestress()[0][1];
		setup.ext_stress_coeff_00 = element->ext_stress_coeff()[0][0];
		setup.ext_stress_coeff_11 = element->ext_stress_coeff()[1][1];
		setup.ext_stress_coeff_01 = element->ext_stress_coeff()[0][1];
		auto C = mepls::utils::tensor::standard_rank4_to_voigt(element->C());
		setup.C_00 = C[0][0];
		setup.C_01 = C[0][1];
		setup.C_02 = C[0][2];
		setup.C_11 = C[1][1];
		setup.C_12 = C[1][2];
		setup.C_22 = C[2][2];
		element_setup.push_back(setup);
	}

	H5::CompType mtype(sizeof(ElementSetupRow));
	//       mtype.insertMember( "type", HOFFSET(ElementSetupRow, type), H5::PredType::NATIVE_UINT);
	mtype.insertMember("prestress_00", HOFFSET(ElementSetupRow, prestress_00),
					   H5::PredType::NATIVE_FLOAT);
	mtype.insertMember("prestress_11", HOFFSET(ElementSetupRow, prestress_11),
					   H5::PredType::NATIVE_FLOAT);
	mtype.insertMember("prestress_01", HOFFSET(ElementSetupRow, prestress_01),
					   H5::PredType::NATIVE_FLOAT);
	mtype.insertMember("ext_stress_coeff_00", HOFFSET(ElementSetupRow, ext_stress_coeff_00),
					   H5::PredType::NATIVE_FLOAT);
	mtype.insertMember("ext_stress_coeff_11", HOFFSET(ElementSetupRow, ext_stress_coeff_11),
					   H5::PredType::NATIVE_FLOAT);
	mtype.insertMember("ext_stress_coeff_01", HOFFSET(ElementSetupRow, ext_stress_coeff_01),
					   H5::PredType::NATIVE_FLOAT);
	mtype.insertMember("C_00", HOFFSET(ElementSetupRow, C_00), H5::PredType::NATIVE_FLOAT);
	mtype.insertMember("C_01", HOFFSET(ElementSetupRow, C_01), H5::PredType::NATIVE_FLOAT);
	mtype.insertMember("C_02", HOFFSET(ElementSetupRow, C_02), H5::PredType::NATIVE_FLOAT);
	mtype.insertMember("C_11", HOFFSET(ElementSetupRow, C_11), H5::PredType::NATIVE_FLOAT);
	mtype.insertMember("C_12", HOFFSET(ElementSetupRow, C_12), H5::PredType::NATIVE_FLOAT);
	mtype.insertMember("C_22", HOFFSET(ElementSetupRow, C_22), H5::PredType::NATIVE_FLOAT);

	hsize_t d[] = {element_setup.size()};
	H5::DataSpace space(1, d);
	H5::DataSet dataset = file.createDataSet(path.c_str(), mtype, space);

	dataset.write(element_setup.data(), mtype);
}


template<int dim>
inline void snapshots(H5::H5File &file,
					  std::string path,
					  const std::vector<mepls::snapshot::Threshold<dim> > &threshold_snapshots,
					  const std::vector<mepls::snapshot::Stress<dim> > &stress_snapshots,
					  const std::vector<mepls::snapshot::DefGrad<dim> > &def_grad_snapshots,
					  const std::vector<mepls::patches::PatchPropertiesSnapshot<dim> > &patch_prop_snapshots,
					  const std::vector<GlobalPropertiesSnapshot<dim>> &global_properties_snapshots)
{

	if(not H5Lexists(file.getId(), path.c_str(), H5P_DEFAULT))
		file.createGroup(path);

	{   /* --------- write threshold snapshots ----------- */

		if(not H5Lexists(file.getId(), (path + "/slip_thresholds").c_str(), H5P_DEFAULT))
			file.createGroup(path + "/slip_thresholds");

		using DataRow = typename mepls::snapshot::Threshold<dim>::DataRow;
		H5::CompType mtype(sizeof(DataRow));
		mtype.insertMember("element", HOFFSET(DataRow, element), H5::PredType::NATIVE_UINT);
		mtype.insertMember("threshold", HOFFSET(DataRow, threshold), H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("angle", HOFFSET(DataRow, angle), H5::PredType::NATIVE_FLOAT);

		unsigned int n = 0;
		for(auto &snapshot : threshold_snapshots)
		{
			hsize_t d[] = {snapshot.data.size()};
			H5::DataSpace space(1, d);
			H5::DataSet dataset = file.createDataSet(
				path + "/slip_thresholds/" + std::to_string(n++), mtype, space);
			dataset.write(snapshot.data.data(), mtype);

			H5::DataSpace att_space(H5S_SCALAR);
			dataset.createAttribute("desired_target", H5::PredType::NATIVE_DOUBLE, att_space).write(
				H5::PredType::NATIVE_DOUBLE, &snapshot.desired_target);
			dataset.createAttribute("recorded_target", H5::PredType::NATIVE_DOUBLE,
									att_space).write(H5::PredType::NATIVE_DOUBLE,
													 &snapshot.recorded_target);
			dataset.createAttribute("index", H5::PredType::NATIVE_UINT, att_space).write(
				H5::PredType::NATIVE_UINT, &snapshot.output_index);
			H5::StrType strdatatype(H5::PredType::C_S1, 32);
			H5std_string strwritebuf_mon(snapshot.monitor_name);
			dataset.createAttribute("monitor_magnitude", strdatatype, att_space).write(strdatatype,
																					   strwritebuf_mon);
			H5std_string strwritebuf_mag(snapshot.recorded_mag);
			dataset.createAttribute("recorded_magnitude", strdatatype, att_space).write(strdatatype,
																						strwritebuf_mag);
		}
	}

	{   /* --------- write stress snapshots ----------- */

		if(not H5Lexists(file.getId(), (path + "/stress").c_str(), H5P_DEFAULT))
			file.createGroup(path + "/stress");

		using DataRow = typename mepls::snapshot::Stress<dim>::DataRow;
		H5::CompType mtype(sizeof(DataRow));
		mtype.insertMember("00", HOFFSET(DataRow, xx), H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("11", HOFFSET(DataRow, yy), H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("01", HOFFSET(DataRow, xy), H5::PredType::NATIVE_FLOAT);

		unsigned int n = 0;
		for(auto &snapshot : stress_snapshots)
		{
			hsize_t d[] = {snapshot.data.size()};
			H5::DataSpace space(1, d);
			H5::DataSet dataset = file.createDataSet(path + "/stress/" + std::to_string(n++), mtype,
													 space);
			dataset.write(snapshot.data.data(), mtype);

			H5::DataSpace att_space(H5S_SCALAR);
			dataset.createAttribute("desired_target", H5::PredType::NATIVE_DOUBLE, att_space).write(
				H5::PredType::NATIVE_DOUBLE, &snapshot.desired_target);
			dataset.createAttribute("recorded_target", H5::PredType::NATIVE_DOUBLE,
									att_space).write(H5::PredType::NATIVE_DOUBLE,
													 &snapshot.recorded_target);
			dataset.createAttribute("index", H5::PredType::NATIVE_UINT, att_space).write(
				H5::PredType::NATIVE_UINT, &snapshot.output_index);
			H5::StrType strdatatype(H5::PredType::C_S1, 32);
			H5std_string strwritebuf_mon(snapshot.monitor_name);
			dataset.createAttribute("monitor_magnitude", strdatatype, att_space).write(strdatatype,
																					   strwritebuf_mon);
			H5std_string strwritebuf_mag(snapshot.recorded_mag);
			dataset.createAttribute("recorded_magnitude", strdatatype, att_space).write(strdatatype,
																						strwritebuf_mag);
		}
	}

	{   /* --------- write deformation gradient snapshots ----------- */

		if(not H5Lexists(file.getId(), (path + "/def_grad").c_str(), H5P_DEFAULT))
			file.createGroup(path + "/def_grad");

		using DataRow = typename mepls::snapshot::DefGrad<dim>::DataRow;
		H5::CompType mtype(sizeof(DataRow));
		mtype.insertMember("00", HOFFSET(DataRow, xx), H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("11", HOFFSET(DataRow, yy), H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("01", HOFFSET(DataRow, xy), H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("10", HOFFSET(DataRow, yx), H5::PredType::NATIVE_FLOAT);

		unsigned int n = 0;
		for(auto &snapshot : def_grad_snapshots)
		{
			hsize_t d[] = {snapshot.data.size()};
			H5::DataSpace space(1, d);
			H5::DataSet dataset = file.createDataSet(path + "/def_grad/" + std::to_string(n++),
													 mtype, space);
			dataset.write(snapshot.data.data(), mtype);

			H5::DataSpace att_space(H5S_SCALAR);
			dataset.createAttribute("desired_target", H5::PredType::NATIVE_DOUBLE, att_space).write(
				H5::PredType::NATIVE_DOUBLE, &snapshot.desired_target);
			dataset.createAttribute("recorded_target", H5::PredType::NATIVE_DOUBLE,
									att_space).write(H5::PredType::NATIVE_DOUBLE,
													 &snapshot.recorded_target);
			dataset.createAttribute("index", H5::PredType::NATIVE_UINT, att_space).write(
				H5::PredType::NATIVE_UINT, &snapshot.output_index);
			H5::StrType strdatatype(H5::PredType::C_S1, 32);
			H5std_string strwritebuf(snapshot.monitor_name);
			dataset.createAttribute("monitor_magnitude", strdatatype, att_space).write(strdatatype,
																					   strwritebuf);
			H5std_string strwritebuf_mag(snapshot.recorded_mag);
			dataset.createAttribute("recorded_magnitude", strdatatype, att_space).write(strdatatype,
																						strwritebuf_mag);
		}
	}

	{   /* --------- write patch snapshots ----------- */

		if(not H5Lexists(file.getId(), (path + "/patches").c_str(), H5P_DEFAULT))
			file.createGroup(path + "/patches");

		using DataRow = typename mepls::patches::PatchPropertiesSnapshot<dim>::DataRow;
		H5::CompType mtype(sizeof(DataRow));
		mtype.insertMember("ref_element", HOFFSET(DataRow, ref_element), H5::PredType::NATIVE_UINT);
		mtype.insertMember("theta", HOFFSET(DataRow, theta), H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("failed", HOFFSET(DataRow, failed), H5::PredType::NATIVE_UINT);
		mtype.insertMember("ss_00", HOFFSET(DataRow, stress_ss_00), H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("ss_11", HOFFSET(DataRow, stress_ss_11), H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("ss_01", HOFFSET(DataRow, stress_ss_01), H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("oi_00", HOFFSET(DataRow, stress_oi_00), H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("oi_11", HOFFSET(DataRow, stress_oi_11), H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("oi_01", HOFFSET(DataRow, stress_oi_01), H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("oi_eps", HOFFSET(DataRow, shear_strain_oi), H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("ee_00", HOFFSET(DataRow, stress_ee_00), H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("ee_11", HOFFSET(DataRow, stress_ee_11), H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("ee_01", HOFFSET(DataRow, stress_ee_01), H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("ss_pe_el", HOFFSET(DataRow, energy_el_ss), H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("ss_pe_conf", HOFFSET(DataRow, energy_conf_ss), H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("x", HOFFSET(DataRow, x), H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("y", HOFFSET(DataRow, y), H5::PredType::NATIVE_FLOAT);


		unsigned int n = 0;
		for(auto &snapshot : patch_prop_snapshots)
		{
			hsize_t d[] = {snapshot.data.size()};
			H5::DataSpace space(1, d);
			H5::DataSet dataset = file
				.createDataSet(path + "/patches/" + std::to_string(n++), mtype, space);
			dataset.write(snapshot.data.data(), mtype);

			H5::DataSpace att_space(H5S_SCALAR);
			dataset.createAttribute("desired_target", H5::PredType::NATIVE_DOUBLE, att_space)
				   .write(H5::PredType::NATIVE_DOUBLE, &snapshot.desired_target);
			dataset.createAttribute("recorded_target", H5::PredType::NATIVE_DOUBLE, att_space)
				   .write(H5::PredType::NATIVE_DOUBLE, &snapshot.recorded_target);
			dataset.createAttribute("index", H5::PredType::NATIVE_UINT, att_space)
				   .write(H5::PredType::NATIVE_UINT, &snapshot.output_index);
			dataset.createAttribute("N_patch", H5::PredType::NATIVE_UINT, att_space)
				   .write(H5::PredType::NATIVE_UINT, &snapshot.N);
			H5::StrType strdatatype(H5::PredType::C_S1, 32);
			H5std_string strwritebuf(snapshot.monitor_name);
			dataset.createAttribute("monitor_magnitude", strdatatype, att_space)
				   .write(strdatatype, strwritebuf);
			H5std_string strwritebuf_mag(snapshot.recorded_mag);
			dataset.createAttribute("recorded_magnitude", strdatatype, att_space)
				   .write(strdatatype, strwritebuf_mag);
		}
	}

	{   /* --------- write full-system properties ----------- */

		if(not H5Lexists(file.getId(), (path + "/global_properties").c_str(), H5P_DEFAULT))
			file.createGroup(path + "/global_properties");

		using DataRow = typename GlobalPropertiesSnapshot<dim>::DataRow;
		H5::CompType mtype(sizeof(DataRow));
		mtype.insertMember("n", HOFFSET(DataRow, n), H5::PredType::NATIVE_UINT);
		mtype.insertMember("C_00", HOFFSET(DataRow, C_00), H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("C_01", HOFFSET(DataRow, C_01), H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("C_02", HOFFSET(DataRow, C_02), H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("C_10", HOFFSET(DataRow, C_10), H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("C_11", HOFFSET(DataRow, C_11), H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("C_12", HOFFSET(DataRow, C_12), H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("C_20", HOFFSET(DataRow, C_20), H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("C_21", HOFFSET(DataRow, C_21), H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("C_22", HOFFSET(DataRow, C_22), H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("theta", HOFFSET(DataRow, theta), H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("oi_eps", HOFFSET(DataRow, oi_eps), H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("ss_00", HOFFSET(DataRow, ss_00), H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("ss_11", HOFFSET(DataRow, ss_11), H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("ss_01", HOFFSET(DataRow, ss_01), H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("oi_00", HOFFSET(DataRow, oi_00), H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("oi_11", HOFFSET(DataRow, oi_11), H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("oi_01", HOFFSET(DataRow, oi_01), H5::PredType::NATIVE_FLOAT);

		unsigned int n = 0;
		for(auto &snapshot : global_properties_snapshots)
		{
			hsize_t d[] = {snapshot.data.size()};
			H5::DataSpace space(1, d);
			H5::DataSet dataset = file.createDataSet(
				path + "/global_properties/" + std::to_string(n++), mtype, space);
			dataset.write(snapshot.data.data(), mtype);

			H5::DataSpace att_space(H5S_SCALAR);
			dataset.createAttribute("desired_target", H5::PredType::NATIVE_DOUBLE, att_space).write(
				H5::PredType::NATIVE_DOUBLE, &snapshot.desired_target);
			dataset.createAttribute("recorded_target", H5::PredType::NATIVE_DOUBLE,
									att_space).write(H5::PredType::NATIVE_DOUBLE,
													 &snapshot.recorded_target);
			dataset.createAttribute("index", H5::PredType::NATIVE_UINT, att_space).write(
				H5::PredType::NATIVE_UINT, &snapshot.output_index);
			H5::StrType strdatatype(H5::PredType::C_S1, 32);
			H5std_string strwritebuf(snapshot.monitor_name);
			dataset.createAttribute("monitor_magnitude", strdatatype, att_space).write(strdatatype,
																					   strwritebuf);
			H5std_string strwritebuf_mag(snapshot.recorded_mag);
			dataset.createAttribute("recorded_magnitude", strdatatype, att_space).write(strdatatype,
																						strwritebuf_mag);
		}
	}

}


inline std::string make_filename(const parameters::Parameters &p)
{
	std::ostringstream file_descriptor_ostrg;
	std::ostringstream L;
	std::ostringstream coupling_constant;
	std::ostringstream elastic_properties;
	std::ostringstream k;
	std::ostringstream nslip;
	std::ostringstream lambda;
	std::ostringstream temperature;
	std::ostringstream beta;
	std::ostringstream gamma;
	std::ostringstream seed;

	L << "L_" << p.sim.Nx << "x" << p.sim.Ny;
	coupling_constant << "+c_" << std::fixed << std::setprecision(6) << p.mat.coupling_constant;
	nslip << "+nslip_" << p.mat.n_slip_systems;
	gamma << "+gamma_" << p.mat.gamma_pl_trans;
	beta << "+beta_" << p.mat.beta;
	elastic_properties << "+G_" << std::fixed << std::setprecision(3) << p.mat.G_quench;
	lambda << "+lambda_" << std::fixed << std::setprecision(4) << p.mat.lambda_quench;
	k << "+k_" << std::fixed << std::setprecision(4) << p.mat.k_quench;
	temperature << "+Tl_" << std::fixed << std::setprecision(2) << p.mat.temperature_liquid;
	seed << "+seed_" << p.sim.seed;

	file_descriptor_ostrg << L.str() << lambda.str() << k.str() << beta.str() << gamma.str()
						  << elastic_properties.str() << coupling_constant.str()
						  << temperature.str() << seed.str();

	return file_descriptor_ostrg.str();
}

} // namespace write



template<int dim>
bool equilibrate_structure_by_rejection(mepls::system::System<dim> &system,
											 const parameters::Parameters &p,
											 unsigned int n_max = 10000)
{
	// this ensures that initially all the thresholds are above the local stress

	if(p.out.verbosity and omp_get_thread_num() == 0)
		std::cout << ">> Equilibrating initial structure... " << std::endl;

	bool relaxed = true;

	for(auto &element : system)
	{
		bool unstable = true;
		unsigned int n = 0;
		while(unstable and n < n_max)
		{
			element->renew_structural_properties();

			unstable = false;
			for(auto &slip : *element)
				if(slip->barrier <= 0)
				{
					unstable = true;
					break;
				}

			++n;
		}

		if(n>=n_max)
		{
			relaxed = false;
			break;
		}
	}


	if(p.out.verbosity and omp_get_thread_num() == 0)
		if(relaxed)
			std::cout << ">>>> Relax. by rejection SUCCESSFUL" << std::endl;
		else
			std::cout << ">>>> Relax. by rejection FAILED" << std::endl;

	return relaxed;
}


template<int dim>
void apply_initial_eigenstrain(mepls::system::System<dim> &system,
								const parameters::Parameters &p)
{
	auto &elements = system.elements;
	auto &solver = system.solver;
	auto &generator = system.generator;

	std::vector<dealii::SymmetricTensor<dim, 2>> eigenstrain(p.sim.Nx * p.sim.Ny);

	std::normal_distribution<double> normal_dist_dev(0., p.mat.init_eigenstrain_std_dev);

	// average volumetric eigenstrain, which is also a random var. since if fluctuates
	// from sample to sample
	std::normal_distribution<double> normal_dist_av_vol(p.mat.init_eigenstrain_av_vol, p.mat.init_eigenstrain_std_av_vol);
	double av_vol_eigenstrain = normal_dist_av_vol(generator);

	std::normal_distribution<double> normal_dist_vol(av_vol_eigenstrain, p.mat.init_eigenstrain_std_vol);

	dealii::SymmetricTensor<dim, 2> eigenstrain_dev_0;
	dealii::SymmetricTensor<dim, 2> eigenstrain_dev_1;
	dealii::SymmetricTensor<dim, 2> eigenstrain_vol;
	dealii::SymmetricTensor<dim, 2> eigenstrain_;

	for(unsigned int n = 0; n < eigenstrain.size(); ++n)
	{
		eigenstrain_dev_0[0][1] = normal_dist_dev(generator);

		eigenstrain_dev_1[0][0] = normal_dist_dev(generator);
		eigenstrain_dev_1[1][1] = -eigenstrain_dev_1[0][0];

		eigenstrain_vol[0][0] = normal_dist_vol(generator);
		eigenstrain_vol[1][1] = eigenstrain_vol[0][0];

		eigenstrain[n] = eigenstrain_dev_0 + eigenstrain_dev_1 + eigenstrain_vol;
	}

	/* ---------------- add eigenstrain and get stress -----------------*/
	solver.clear();

	for(unsigned int n = 0; n < eigenstrain.size(); ++n)
		solver.add_eigenstrain(n, eigenstrain[n]);

	// dummy event to make the system update the elastic state of the elements
	mepls::event::Driving<dim> init_eigenstrain_event;
	init_eigenstrain_event.activation_protocol = mepls::dynamics::Protocol::prestress;
	system.add(init_eigenstrain_event);

	// renew the properties, so they take into account the stress.
	// This matters if there is pressure sensitivity. This call ensures
	// that all the elements are renewed at least once, since the structure
	// equilibration done after this, will renew only the unstable ones
	for(auto &element : elements)
		element->renew_structural_properties();

	// use rejection instead of relaxation, because we don't want the
	// stress field to be altered anymore. Note: depending on the stress
	// field and the threshold distribution, this method might never find
	// a stable configuration. In that case, relaxation is mandatory.
	bool relaxed = equilibrate_structure_by_rejection(system, p, 1000000);

	if(not relaxed)
	{
		mepls::utils::ContinueSimulation continue_simulation;
		mepls::dynamics::relaxation(system, continue_simulation);
		if(not continue_simulation())
		{
			std::cout << continue_simulation << std::endl;
			abort();
		}
	}
}


template<int dim>
void simulate_parent_liquid_KMC(mepls::system::System<dim> &system,
							mepls::history::History<dim> &history,
							const parameters::Parameters &p,
							mepls::utils::ContinueSimulation &continue_simulation)
{
	mepls::dynamics::KMC<dim> kmc;
	auto &macro_evolution = history.macro_evolution;

	std::vector<double> rolling_av_stress;
	double rolling_av_stress_old = 0.;
	double rolling_av_stress_new = 0.;
	mepls::utils::ContinueSimulation continue_kmc;
	unsigned int i = 0;
	unsigned int n = 100;

	while(continue_kmc())
	{
		++i;

		kmc(system);
		mepls::dynamics::relaxation(system, continue_kmc);

		if(p.out.verbosity and omp_get_thread_num() == 0)
			std::cout << i << std::endl;

		history.add_macro(system);

		if(i % n == 0 and i > 1000)
		{
			rolling_av_stress_new = 0.;
			for(unsigned int j = macro_evolution.size() - n; j < macro_evolution.size();
				++j)
				rolling_av_stress_new += macro_evolution[j].av_vm_stress;
			rolling_av_stress_new /= double(n);

			continue_kmc( std::abs(rolling_av_stress_new-rolling_av_stress_old)
			/rolling_av_stress_old > 0.01, "KMC reached the stationary state" );

			rolling_av_stress_old = rolling_av_stress_new;
		}
	}
}


template<int dim>
void simulate_parent_liquid_MH(mepls::system::System<dim> &system,
							mepls::history::History<dim> &history,
							const parameters::Parameters &p,
							mepls::utils::ContinueSimulation &continue_simulation)
{
	mepls::dynamics::MetropolisHastings<dim> mh;
	mh.T = p.mat.temperature_liquid;

	auto &macro_evolution = history.macro_evolution;

	std::vector<double> rolling_av_energy;
	double rolling_av_energy_old = 0.;
	double rolling_av_energy_new = 0.;
	mepls::utils::ContinueSimulation continue_mh;
	unsigned int i_accepted = 0;
	unsigned int i_total = 0;
	unsigned int n = 200;

	while(i_accepted < 2000)//continue_mh())
	{
		bool accepted = mh(system, true);

		if(accepted)
		{
			++i_accepted;

			history.add_macro(system);

			if(p.out.verbosity and omp_get_thread_num() == 0)
				std::cout << i_accepted << " " << system.macrostate["ext_stress"] << std::endl;

			if(i_accepted % n == 0 and i_accepted > 1000)
			{

				rolling_av_energy_new = 0.;
				for(unsigned int j = macro_evolution.size() - n; j < macro_evolution.size();
					++j)
					rolling_av_energy_new += macro_evolution[j].av_energy_conf +
											 macro_evolution[j].av_energy_el;
				rolling_av_energy_new /= double(n);

				if(rolling_av_energy_old!=0.)
				{
					continue_mh( std::abs(rolling_av_energy_new/rolling_av_energy_old-1.) > 5e-3, "MH reached the stationary state" );

					rolling_av_energy_old = rolling_av_energy_new;
				}else{
					// initialization of rolling_av_energy_old
					rolling_av_energy_old = rolling_av_energy_new;
				}

			}
		}

		++i_total;
	}
}


template<int dim>
std::vector<element::Anisotropic<dim> *> create_elements(const parameters::Parameters &p,
														 std::mt19937 &generator)
{
	std::vector<element::Anisotropic<dim> *> elements;

	for(double n = 0; n < p.sim.Nx * p.sim.Ny; ++n)
	{
		typename element::Anisotropic<dim>::Config conf;

		conf.G = p.mat.G;
		conf.K = p.mat.K;
		conf.G_quench = p.mat.G_quench;
		conf.K_quench = p.mat.K_quench;
		conf.alpha_tau = p.mat.alpha_tau;
		conf.k_quench = p.mat.k_quench;
		conf.k = p.mat.k;
		conf.lambda_quench = p.mat.lambda_quench;
		conf.lambda = p.mat.lambda;
		conf.n_slip_systems = p.mat.n_slip_systems;
		conf.coupling_constant = p.mat.coupling_constant;
		conf.activation_rate_0 = p.mat.activation_rate;
		conf.number = n;
		conf.gamma_pl_trans = p.mat.gamma_pl_trans;
		conf.beta = p.mat.beta;
		conf.A = p.mat.A;
		conf.B = p.mat.B;
		conf.av_U = p.mat.av_U;
		conf.std_U = p.mat.std_U;
		conf.temperature = p.mat.temperature_liquid;

		elements.push_back(
			new element::Anisotropic<dim>(generator, conf));
	}

	return elements;
}


template<int dim>
void perform_reloading(mepls::system::System<dim> &system,
					   mepls::history::History<dim> &history,
					   bool is_forward,
					   const parameters::Parameters &p)
{
	// Simulate using a copy of the original system; Clear the solver, and use the element stress
	// as prestress. Clean the deformation of the elements;
	system.solver.clear();

	mepls::element::Vector<dim> elements_replica;
	for(auto &element : system)
	{
		auto element_copy = element->make_copy();
		element_copy->state_to_prestress();
		elements_replica.push_back(element_copy);
	}

	auto system_replica = system.get_new_instance(elements_replica, system.solver,
												  system.generator);
	auto &macrostate = system_replica->macrostate;
	system_replica->set_history(history);

	// initiate dynamics using copied system
	mepls::utils::ContinueSimulation continue_loading;
	history.add_macro( *system_replica );

	while(continue_loading())
	{
		if(p.out.verbosity and omp_get_thread_num() == 0)
			std::cout << history.index() << " | " << std::fixed << macrostate["total_strain"]
					  << " " << macrostate["ext_stress"] << " " << macrostate["pressure"]
					  << std::endl;

		mepls::dynamics::finite_extremal_dynamics_step(1e-4 * 0.5, *system_replica, is_forward);
		history.add_macro( *system_replica );

		mepls::dynamics::relaxation(*system_replica, continue_loading);
		history.add_macro( *system_replica );

		continue_loading(std::abs(macrostate["total_strain"]) < 0.4 / 2., "total_strain limit reached");
	}

	delete system_replica;

	if(p.out.verbosity and omp_get_thread_num() == 0)
	{
		std::cout << continue_loading << std::endl;
		std::cout << "Reloading " << (is_forward ? "forward" : "backward") << " finished"
				  << std::endl;
	}
}


template<int dim>
double elastic_properties_evolution(mepls::system::System<dim> &system, 
									mepls::elasticity_solver::LeesEdwards<dim> &solver, 
									double G_old, 
									double G_stat,
									mepls::utils::ContinueSimulation &continue_simulation)
{
	auto &elements = system.elements;

	// reassemble the elastic properties if the change in the shear
	// modulus is big enough and if it's different enough from the stationary
	// value
	dealii::SymmetricTensor<4,dim> av_C;
	for(auto &element : elements)
		av_C += element->C();
	av_C /= double(elements.size());

	double G_new = av_C[0][1][0][1];

	if(std::abs(G_new/G_old-1)>0.001 and std::abs(G_old/G_stat-1)>0.005)
	{
		G_old = G_new;

		solver.reassemble_with_new_stiffness(av_C);

		// we call add which will inform the macrostate and the elements
		// about the change in external stress due to the change in global
		// stiffness and will also record that change in the history
		const std::vector<mepls::event::Plastic<dim>> added_yielding;
		mepls::event::Driving<dim> driving_event_variation_stiffness;
		driving_event_variation_stiffness.activation_protocol = mepls::dynamics::Protocol::variation_stiffness;
		system.add(driving_event_variation_stiffness);

		auto state = solver.get_state();
		mepls::element::calculate_local_stress_coefficients_central(elements, solver);
		mepls::element::calculate_ext_stress_coefficients(elements, solver);
		solver.set_state(state);

		// relax to ensure there are no unstable elements after the change
		// in the elastic properties (if the shear modulus rises with strain,
		// that will lead to stress rises that can unstabilise elements)
		mepls::dynamics::relaxation(system, continue_simulation);
	}

	return G_old;
}

} // namespace espci


#endif //_ESPCI_H

