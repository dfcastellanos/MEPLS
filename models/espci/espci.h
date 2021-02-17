// Copyright 2020 David Fernández Castellanos
// License: BSD-3
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
//    this list of conditions and the following disclaimer.
//
//2. Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
//
//3. Neither the name of the copyright holder nor the names of its contributors
//   may be used to endorse or promote products derived from this software
//   without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

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
	// prestress
	bool prestress = true;
	double prestress_std_shear = 0.117;
	double prestress_std_pressure = 0.1;
	double prestress_av_pressure = 0.1;
	double prestress_std_av_pressure = 0.2;

	// threshold distribution
	double lambda = 5.;
	double lambda_quench = 5.;
	double k = 1.9;
	double k_quench = 1.9;
	double alpha_tau = 0.;

	// modulus distribution
	double average_G = 1.;
	double average_G_quench = 1.;
	double average_K = 1.;
	double weibull_shape_G = 1.;
	double weibull_shape_K = 1.;

	// others
	double coupling_constant = 0.1;
	double temperature = 1e-2;
	double activation_rate = 1.;
	unsigned int n_slip_systems = 1;
	unsigned int n_slip_systems_quench = 1;

	void declare_entries(dealii::ParameterHandler &prm)
	{
		/*! Declare the variables of this struct into the input parser object. */

		prm.enter_subsection("Material");

		prm.declare_entry("prestress_std_shear", mepls::utils::str::to_string(prestress_std_shear),
						  dealii::Patterns::Double(), "");
		prm.declare_entry("prestress_std_pressure",
						  mepls::utils::str::to_string(prestress_std_pressure),
						  dealii::Patterns::Double(), "");
		prm.declare_entry("prestress_av_pressure",
						  mepls::utils::str::to_string(prestress_av_pressure),
						  dealii::Patterns::Double(), "");
		prm.declare_entry("prestress_std_av_pressure",
						  mepls::utils::str::to_string(prestress_std_av_pressure),
						  dealii::Patterns::Double(), "");

		prm.declare_entry("average_G", mepls::utils::str::to_string(average_G),
						  dealii::Patterns::Double(), "");
		prm.declare_entry("average_G_quench", mepls::utils::str::to_string(average_G_quench),
						  dealii::Patterns::Double(), "");
		prm.declare_entry("average_K", mepls::utils::str::to_string(average_K),
						  dealii::Patterns::Double(), "");
		prm.declare_entry("weibull_shape_G", mepls::utils::str::to_string(weibull_shape_G),
						  dealii::Patterns::Double(), "");
		prm.declare_entry("weibull_shape_K", mepls::utils::str::to_string(weibull_shape_K),
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
		prm.declare_entry("temperature", mepls::utils::str::to_string(temperature),
						  dealii::Patterns::Double(0.0), "");
		prm.declare_entry("activation_rate", mepls::utils::str::to_string(activation_rate),
						  dealii::Patterns::Double(0.0), "");

		prm.declare_entry("n_slip_systems", mepls::utils::str::to_string(n_slip_systems),
						  dealii::Patterns::Integer(0), "");
		prm.declare_entry("n_slip_systems_quench",
						  mepls::utils::str::to_string(n_slip_systems_quench),
						  dealii::Patterns::Integer(0), "");
		prm.declare_entry("prestress", mepls::utils::str::to_string(prestress),
						  dealii::Patterns::Bool(), "");

		prm.leave_subsection();

	}

	void load(dealii::ParameterHandler &prm)
	{
		/*! Load values parsed by the input parser object into the structs variables. */

		prm.enter_subsection("Material");

		prestress_std_shear = prm.get_double("prestress_std_shear");
		prestress_std_pressure = prm.get_double("prestress_std_pressure");
		prestress_av_pressure = prm.get_double("prestress_av_pressure");
		prestress_std_av_pressure = prm.get_double("prestress_std_av_pressure");

		average_G = prm.get_double("average_G");
		average_G_quench = prm.get_double("average_G_quench");
		average_K = prm.get_double("average_K");
		weibull_shape_G = prm.get_double("weibull_shape_G");
		weibull_shape_K = prm.get_double("weibull_shape_K");

		lambda = prm.get_double("lambda");
		lambda_quench = prm.get_double("lambda_quench");
		k = prm.get_double("k");
		k_quench = prm.get_double("k_quench");
		alpha_tau = prm.get_double("alpha_tau");
		coupling_constant = prm.get_double("coupling_constant");
		temperature = prm.get_double("temperature");
		activation_rate = prm.get_double("activation_rate");

		n_slip_systems = prm.get_integer("n_slip_systems");
		n_slip_systems_quench = prm.get_integer("n_slip_systems_quench");
		prestress = prm.get_bool("prestress");

		prm.leave_subsection();
	}

};


/*! Struct with parameters related to the simulation setup. */
struct Simulation
{
	unsigned int n_theta = 1;
	std::vector<int> N_probe_list;
	unsigned int n_rep = 1;

	unsigned int Nx = 32;
	/*!<  Number of elements in the x-direction of the rectangular lattice forming the system. */

	unsigned int Ny = 32;
	/*!< Number of elements in the y-direction of the rectangular lattice forming the system.  */

	unsigned int seed = 1234567;
	/*!<  Seed to initialized the random number engine. */

	double initial_load = 0.;
	/*!< Value of the load applied at the beginning of the simulation.  */

	std::string loading_mode = "pure_shear";
	/*!< Loading conditions, e.g., pure shear, compression etc. */

	std::string control_mode = "displacement";
	/*!< Mode of driving the system. It can be displacement or traction controlled (equivalent to
	 * strain or stress controlled). */

	std::string trigger = "extremal_dynamics";
	/*!<  Protocol (see \ref dynamics) used to unstabilize the system in the main simulation loop.
	 * Posibilities are kmc, extremal_dynamics or finite_extremal_dynamics. */

	std::string monitor_name = "total_strain";
	/*!< Magnitude used to check whether the simulation should stop and whether a snapshot should be
	 *  taken (see \ref snapshot). The posibilities correspond to the keys in
	 *  \ref MacroState.monitor_map. */

	double monitor_limit = 3.0;
	/*!<  Value of the magnitude defined by \ref monitor_name at which the simulation must stop. */

	double fracture_limit = 2;
	/*!<  Number of plastic events triggered during a call to \ref mepls::dynamics::relaxation
	 * divided by the number of elements composing the system at which the simulation must stop.
	 * A default value of 2 means that during a relaxation, every element has deformed twice,
	 * which is very big. Such value can, in most cases, indicate that the material has fractured
	 *  and therefore the simulation must stop. */

	double std_blur = 0.;
	/*!< Value of the standard deviation (in units of the mesoscale element linear length) of the
	 * Gaussian kernel used for blurring the strain field calculated from the elasticity solver. */

	bool kmc_quench = true;
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
		prm.declare_entry("N_probe_list", "", dealii::Patterns::FileName(), "");
		prm.declare_entry("Nx", mepls::utils::str::to_string(Nx), dealii::Patterns::Integer(0), "");
		prm.declare_entry("Ny", mepls::utils::str::to_string(Ny), dealii::Patterns::Integer(0), "");
		prm.declare_entry("seed", mepls::utils::str::to_string(seed), dealii::Patterns::Integer(0),
						  "");
		prm.declare_entry("loading_mode", loading_mode, dealii::Patterns::Selection(
			"pure_shear|pure_shear_pbc|simple_shear|compression"), "");
		prm.declare_entry("initial_load", mepls::utils::str::to_string(initial_load),
						  dealii::Patterns::Double(0.0), "");
		prm.declare_entry("control_mode", control_mode,
						  dealii::Patterns::Selection("displacement|traction"), "");
		prm.declare_entry("trigger", trigger, dealii::Patterns::Selection(
			"extremal_dynamics|kmc|finite_extremal_dynamics"), "");
		prm.declare_entry("monitor_name", monitor_name, dealii::Patterns::Selection(
			"av_plastic_strain|load|ext_stress|time|total_strain"), "");
		prm.declare_entry("monitor_limit", mepls::utils::str::to_string(monitor_limit),
						  dealii::Patterns::Double(0.0), "");
		prm.declare_entry("fracture_limit", mepls::utils::str::to_string(fracture_limit),
						  dealii::Patterns::Double(0.0), "");
		prm.declare_entry("std_blur", mepls::utils::str::to_string(std_blur),
						  dealii::Patterns::Double(0.0), "");
		prm.declare_entry("kmc_quench", mepls::utils::str::to_string(kmc_quench),
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
		loading_mode = prm.get("loading_mode");
		initial_load = prm.get_double("initial_load");
		control_mode = prm.get("control_mode");
		trigger = prm.get("trigger");
		monitor_name = prm.get("monitor_name");
		monitor_limit = prm.get_double("monitor_limit");
		fracture_limit = prm.get_double("fracture_limit");
		std_blur = prm.get_double("std_blur");
		N_probe_list = mepls::utils::str::parse_list_integers(prm.get("N_probe_list"));
		kmc_quench = prm.get_bool("kmc_quench");
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
	 *  threshold, stress, def_grad, local_probe. */

	double snapshots_min = 0.;
	/*!< Minimum value of the magnitude \ref Simulation.monitor_name from which snapshots should
	 * start to be taken (see \ref snapshot).  */

	double snapshots_max = 1.;
	/*!<  Maximum value of the magnitude \ref Simulation.monitor_name until which snapshots should
	 * be taken (see \ref snapshot). */

	double snapshots_interval = 0.05;
	/*!< Interval to create a list of value from \ref snapshots_min to \ref snapshots_max at which
	 *  snapshots should be taken (see \ref snapshot).   */

	double snapshots_sensitivity = 0.01;

	/*!< Numerical tolerance to decide whether the value of \ref Simulation.monitor_name is that at
	 * which a snapshot must be taken (see \ref snapshot). */

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
struct Standard
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
		double coupling_constant = 0.;
		double alpha_tau = 0.;
		double temperature = 1e-2;
		double activation_rate = 1.;
		unsigned int thermal = 0;
	};

	Oriented(std::mt19937 &generator_, const Config &conf_)
		:
		mepls::slip::Slip<dim>(),
		conf(conf_),
		threshold_0(conf_.threshold),
		unif_distribution(0, 1.),
		normal_dist(0, 1),
		generator(generator_)
	{
		M_Assert(threshold_0 > 0., "Expected threshold > 0");

		angle = mepls::utils::mod(conf.angle, M_PI);
		M = mepls::utils::tensor::make_schmid<dim>(angle);M_Assert(M[0][1] == M[1][0], "");M_Assert(
			M[0][0] == -M[1][1], "");
	};

	void update() override
	{
		eff_shear_stress = M * parent->stress();
		pressure = dealii::trace(parent->stress()) / double(dim);
		add_pressure_effects();
		modify_barrier();

		M_Assert(threshold > 0., "Expected threshold > 0");M_Assert(not std::isnan(threshold),
																	"Threshold is NaN");M_Assert(
			not std::isnan(eff_shear_stress), "eff_shear_stress is NaN");
	}

	void add_pressure_effects()
	{
		threshold = threshold_0 - conf.alpha_tau * pressure;
		if(threshold < 0.)
			threshold = 1e-3;
	}

	void modify_barrier()
	{
		barrier = threshold - eff_shear_stress;
		activation_rate = conf.activation_rate * std::exp(-barrier / conf.temperature);
	}

	double get_critical_load_increment() override
	{
		M_Assert(barrier >= 0,
				 "Negative barrier found while looking for critical load increment.");M_Assert(
			(parent->ext_stress_coeff()[0][0] != 0.) or (parent->ext_stress_coeff()[1][1] != 0.) or (parent->ext_stress_coeff()[0][1] != 0.),
			"Ext. stress coeffs. are all 0");

		return barrier / (M * parent->ext_stress_coeff() + 0.5 * conf.alpha_tau * dealii::trace(
			parent->ext_stress_coeff()));
	}

	dealii::SymmetricTensor<2, dim> get_eigenstrain_increment() override
	{
		std::gamma_distribution<double> g_a(1, 1);
		std::gamma_distribution<double> g_b(1 - conf.coupling_constant, 1);
		double za = g_a(generator);
		double zb = g_b(generator);
		double eff_shear_stress_variation = -(eff_shear_stress + conf.thermal * (barrier + conf.temperature)) * za / (za + zb);
//		double eff_shear_stress_variation = - (eff_shear_stress + conf.thermal*barrier ) * za/(za+zb);

		M_Assert(eff_shear_stress *eff_shear_stress_variation<1,
					 "Eff. shear stress variation does not reduce local stress. Is the orientation of plastic event according to local stress?");

		double A = M * parent->S() * M;
		double gamma = 2 * eff_shear_stress_variation / A;
		//			double gamma = conf.coupling_constant;
		dealii::SymmetricTensor<2, dim> eigenstrain_dev = gamma * M;

		return eigenstrain_dev;
	}

	mepls::slip::Slip<dim> *make_copy() override
	{
		return new slip::Oriented<dim>(generator, conf);
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
	mutable std::uniform_real_distribution<double> unif_distribution;
	mutable std::normal_distribution<double> normal_dist;
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
		double average_G = 1.;
		double average_K = 1.;
		double weibull_shape_G = 1.;
		double weibull_shape_K = 1.;

		double alpha_tau = 0.;
		unsigned int thermal = 0;

		double coupling_constant = 0.2;
		double temperature = 1e-2;
		double activation_rate = 1.;

		unsigned int n_slip_systems = 1;
		unsigned int number = 0;
		dealii::SymmetricTensor<2, dim> prestress;
	};


	Anisotropic(std::piecewise_linear_distribution<double> &threshold_distribution_,
				std::mt19937 &generator_,
				const Config &conf_)
		:
		mepls::element::Element<dim>(),
		threshold_distribution_ptr(&threshold_distribution_),
		generator(generator_),
		unif_distribution(0, 1),
		conf(conf_)
	{
		this->number(conf.number);
		this->prestress(conf.prestress);

		renew_elastic_properties();
		renew_thresholds();
	}


	Anisotropic(Anisotropic *input_element)
		:
		mepls::element::Element<dim>(),
		threshold_distribution_ptr(input_element->threshold_distribution_ptr),
		generator(input_element->generator),
		unif_distribution(0, 1),
		conf(input_element->conf)
	{
	}


	void renew_structural_properties_impl(mepls::element::RenewInstruct<dim> &renew_instruct) override
	{
		if(renew_instruct.elastic_properties)
			renew_elastic_properties();

		if(renew_instruct.slip_properties)
			renew_thresholds();
	}


	void renew_elastic_properties()
	{
		double G = mepls::utils::rand::get_weibull_rand_av(conf.weibull_shape_G, conf.average_G,
														   unif_distribution(generator));
		auto C = mepls::utils::tensor::make_mandel_anisotropic_stiffness<dim>(conf.average_K, G, G,
																			  0.);
		this->C(mepls::utils::tensor::mandel_to_standard_rank4<dim>(C));
	}


	void renew_thresholds()
	{
		this->remove_slip_systems();

		typename slip::Oriented<dim>::Config slip_conf;
		slip_conf.alpha_tau = conf.alpha_tau;
		slip_conf.coupling_constant = conf.coupling_constant;
		slip_conf.temperature = conf.temperature;
		slip_conf.activation_rate = conf.activation_rate;
		slip_conf.thermal = conf.thermal;

		auto &threshold_distribution = *threshold_distribution_ptr;

		//			slip_conf.angle = 0;
		//			slip_conf.threshold = threshold_distribution(generator);
		//			this->add_slip_system( new slip::Oriented<dim>(generator, slip_conf) );
		//
		//			slip_conf.angle = 0 + M_PI/2.;
		//			this->add_slip_system( new slip::Oriented<dim>(generator, slip_conf) );

		for(unsigned int i = 0; i < conf.n_slip_systems; ++i)
		{
			double alpha2 = unif_distribution(generator) * M_PI;

			slip_conf.angle = alpha2;
			slip_conf.threshold = threshold_distribution(generator);
			this->add_slip_system(new slip::Oriented<dim>(generator, slip_conf));

			slip_conf.angle = alpha2 + M_PI / 2.;
			slip_conf.threshold = threshold_distribution(generator);
			this->add_slip_system(new slip::Oriented<dim>(generator, slip_conf));

			slip_conf.angle = alpha2 + M_PI / 4.;
			slip_conf.threshold = threshold_distribution(generator);
			this->add_slip_system(new slip::Oriented<dim>(generator, slip_conf));

			slip_conf.angle = alpha2 + M_PI / 4. + M_PI / 2.;
			slip_conf.threshold = threshold_distribution(generator);
			this->add_slip_system(new slip::Oriented<dim>(generator, slip_conf));
		}
	}


	mepls::element::Element<dim> *make_copy_impl() override
	{
		return new Anisotropic<dim>(this);
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
			slip_conf.thermal = conf.thermal;
			slip_conf.alpha_tau = conf.alpha_tau;
			slip_conf.coupling_constant = conf.coupling_constant;
			slip_conf.temperature = conf.temperature;
			slip_conf.activation_rate = conf.activation_rate;
			slip_conf.thermal = conf.thermal;
		}
	}

	void threshold_distribution(std::piecewise_linear_distribution<double> *threshold_distribution_ptr_)
	{
		threshold_distribution_ptr = threshold_distribution_ptr_;
	}

protected:
	std::piecewise_linear_distribution<double> *threshold_distribution_ptr;
	std::mt19937 &generator;
	std::uniform_real_distribution<double> unif_distribution;
	Config conf;
};

} // namespace element


namespace history
{


struct MacroSummaryRow
{
	/*! Struct to store system-scale properties */

	float av_vm_plastic_strain = 0.;
	float std_vm_plastic_strain = 0.;
	float av_vm_stress = 0.;
	float av_pressure = 0.;
	float std_vm_stress = 0.;
	float std_pressure = 0.;
	float av_slip_threshold = 0.;
	float std_slip_threshold = 0.;
	double time = 0.;
	float total_strain = 0.;
	float ext_stress = 0.;
	float av_potential_energy = 0.;
	float std_potential_energy = 0.;
};



template<int dim>
class EventAndMacro: public mepls::history::History<dim>
{
  public:

	EventAndMacro(const std::string &name_ = "history")
		:
		mepls::history::History<dim>(),
		name(name_)
	{
		/*! Constructor. */
	}


	void clear()
	{
		macro_evolution.clear();
		mepls::history::History<dim>::clear();
	}


   void add_macro(const mepls::system::System<dim> &system)
   {
      if(closed)
         return;

		MacroSummaryRow data;

		const auto &elements = system.elements;
		const auto &macrostate = system.macrostate;

		double av_vm_plastic_strain2 = 0.;
		double av_vm_stress2 = 0.;
		double av_pressure2 = 0.;
		double av_slip_threshold2 = 0.;
		double av_potential_energy2 = 0.;
		unsigned int n_total_slip = 0;

		for(auto &element : elements)
		{
			double vm_plastic_strain = element->integrated_vm_eigenstrain();
			double vm_stress = mepls::utils::get_von_mises_equivalent_stress(element->stress());
			double pressure = -dealii::trace(element->stress()) / double(dim);
			double potential_energy = 0.5 * dealii::invert(
				element->C()) * element->stress() * element->stress();
			data.av_vm_plastic_strain += vm_plastic_strain;
			data.av_vm_stress += vm_stress;
			data.av_pressure += pressure;
			data.av_potential_energy += potential_energy;
			av_vm_plastic_strain2 += vm_plastic_strain * vm_plastic_strain;
			av_vm_stress2 += vm_stress * vm_stress;
			av_pressure2 += pressure * pressure;
			av_potential_energy2 += potential_energy * potential_energy;
			for(auto &slip : *element)
			{
				data.av_slip_threshold += slip->threshold;
				av_slip_threshold2 += slip->threshold * slip->threshold;
				++n_total_slip;
			}

		}

		double N = double(elements.size());
		data.av_vm_plastic_strain /= N;
		data.av_vm_stress /= N;
		data.av_pressure /= N;
		data.av_potential_energy /= N;
		data.av_slip_threshold /= double(n_total_slip);
		av_vm_plastic_strain2 /= N;
		av_vm_stress2 /= N;
		av_pressure2 /= N;
		av_potential_energy2 /= N;
		av_slip_threshold2 /= double(n_total_slip);
		data.std_vm_plastic_strain = std::sqrt(
			av_vm_plastic_strain2 - data.av_vm_plastic_strain * data.av_vm_plastic_strain);
		data.std_vm_stress = std::sqrt(av_vm_stress2 - data.av_vm_stress * data.av_vm_stress);
		data.std_pressure = std::sqrt(av_pressure2 - data.av_pressure * data.av_pressure);
		data.std_potential_energy = std::sqrt(
			av_potential_energy2 - data.av_potential_energy * data.av_potential_energy);
		data.std_slip_threshold = std::sqrt(
			av_slip_threshold2 - data.av_slip_threshold * data.av_slip_threshold);

		data.time = macrostate["time"];
		data.total_strain = macrostate["total_strain"];
		data.ext_stress = macrostate["ext_stress"];

      macro_evolution.push_back(data);
   }

	std::vector<MacroSummaryRow> macro_evolution;
	std::string name;

	using mepls::history::History<dim>::closed;
	using mepls::history::History<dim>::index;
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
		mepls::patches::PatchPropertiesTensorial<dim> probed_stress;

		// external loading is shear with only epsxy!=0.
		double theta = 0.;

		// to compare with the MD full system global_properties, the external strain discrete
		// incremement is here 1e-4 (different from the 1e-3 used for the patches)
		double dgamma = 1e-4;
		mepls::patches::apply_patch_shear_test<dim>(probed_stress, *system_replica,
													continue_shear_test, false, dgamma);

		global_properties.oi_eps = probed_stress.resolved_elastic_shear_strain_oi;
		global_properties.ss_00 = probed_stress.stress_ss[0][0];
		global_properties.ss_11 = probed_stress.stress_ss[1][1];
		global_properties.ss_01 = probed_stress.stress_ss[0][1];
		global_properties.oi_00 = probed_stress.stress_oi[0][0];
		global_properties.oi_11 = probed_stress.stress_oi[1][1];
		global_properties.oi_01 = probed_stress.stress_oi[0][1];

		solver.clear();


		data.push_back(global_properties);
	};

	std::vector<DataRow> data;
	/*!< Container to store the recorde data. */

	std::string recorded_mag;
	/*!< Name of the field storaged in \ref data. */

	std::string monitor_name;
	/*!< Name of the magnitude used to check whether the snapshot should be taken or not. */

	double desired_target;
	/*!< Value of the \ref monitor_name at which we desired to take the snapshot. */

	double recorded_target;
	/*!< Value of the \ref monitor_name at which the snapshot is actually taken. */

	unsigned int output_index;
	/*!< Global event index from the \ref mepls::event::History at which the snapshot is taken. */
};


namespace write
{

inline void file_attrs(H5::H5File &file, const parameters::Standard &p)
{
	H5::DataSpace att_space(H5S_SCALAR);
	file.createAttribute("lambda", H5::PredType::NATIVE_DOUBLE, att_space).write(
		H5::PredType::NATIVE_DOUBLE, &p.mat.lambda);
	file.createAttribute("k", H5::PredType::NATIVE_DOUBLE, att_space).write(
		H5::PredType::NATIVE_DOUBLE, &p.mat.k);
	file.createAttribute("alpha_tau", H5::PredType::NATIVE_DOUBLE, att_space).write(
		H5::PredType::NATIVE_DOUBLE, &p.mat.alpha_tau);
	file.createAttribute("average_G", H5::PredType::NATIVE_DOUBLE, att_space).write(
		H5::PredType::NATIVE_DOUBLE, &p.mat.average_G);
	file.createAttribute("weibull_shape_G", H5::PredType::NATIVE_DOUBLE, att_space).write(
		H5::PredType::NATIVE_DOUBLE, &p.mat.weibull_shape_G);
	file.createAttribute("average_K", H5::PredType::NATIVE_DOUBLE, att_space).write(
		H5::PredType::NATIVE_DOUBLE, &p.mat.average_K);
	file.createAttribute("weibull_shape_K", H5::PredType::NATIVE_DOUBLE, att_space).write(
		H5::PredType::NATIVE_DOUBLE, &p.mat.weibull_shape_K);
	file.createAttribute("Nx", H5::PredType::NATIVE_UINT, att_space).write(
		H5::PredType::NATIVE_UINT, &p.sim.Nx);
	file.createAttribute("Ny", H5::PredType::NATIVE_UINT, att_space).write(
		H5::PredType::NATIVE_UINT, &p.sim.Ny);
	file.createAttribute("seed", H5::PredType::NATIVE_UINT, att_space).write(
		H5::PredType::NATIVE_UINT, &p.sim.seed);
	file.createAttribute("coupling_constant", H5::PredType::NATIVE_DOUBLE, att_space).write(
		H5::PredType::NATIVE_DOUBLE, &p.mat.coupling_constant);
	file.createAttribute("temperature", H5::PredType::NATIVE_DOUBLE, att_space).write(
		H5::PredType::NATIVE_DOUBLE, &p.mat.temperature);
	file.createAttribute("activation_rate", H5::PredType::NATIVE_DOUBLE, att_space).write(
		H5::PredType::NATIVE_DOUBLE, &p.mat.activation_rate);
	file.createAttribute("prestress_av_pressure", H5::PredType::NATIVE_DOUBLE, att_space).write(
		H5::PredType::NATIVE_DOUBLE, &p.mat.prestress_av_pressure);
	file.createAttribute("prestress_std_av_pressure", H5::PredType::NATIVE_DOUBLE, att_space).write(
		H5::PredType::NATIVE_DOUBLE, &p.mat.prestress_std_av_pressure);
	file.createAttribute("prestress_std_pressure", H5::PredType::NATIVE_DOUBLE, att_space).write(
		H5::PredType::NATIVE_DOUBLE, &p.mat.prestress_std_pressure);
	file.createAttribute("prestress_std_shear", H5::PredType::NATIVE_DOUBLE, att_space).write(
		H5::PredType::NATIVE_DOUBLE, &p.mat.prestress_std_shear);
	file.createAttribute("n_slip_systems", H5::PredType::NATIVE_UINT, att_space).write(
		H5::PredType::NATIVE_UINT, &p.mat.n_slip_systems);

	{
		H5std_string strwritebuf;
		strwritebuf = p.sim.monitor_name;
		H5::StrType strdatatype(H5::PredType::C_S1, strwritebuf.size());
		file.createAttribute("monitor_name", strdatatype, att_space).write(strdatatype,
																		   strwritebuf);
	}
	{
		H5std_string strwritebuf;
		strwritebuf = p.sim.trigger;
		H5::StrType strdatatype(H5::PredType::C_S1, strwritebuf.size());
		file.createAttribute("trigger", strdatatype, att_space).write(strdatatype, strwritebuf);
	}
	{
		H5std_string strwritebuf;
		strwritebuf = p.sim.loading_mode;
		H5::StrType strdatatype(H5::PredType::C_S1, strwritebuf.size());
		file.createAttribute("loading_mode", strdatatype, att_space).write(strdatatype,
																		   strwritebuf);
	}
	{
		H5std_string strwritebuf;
		strwritebuf = p.sim.control_mode;
		H5::StrType strdatatype(H5::PredType::C_S1, strwritebuf.size());
		file.createAttribute("control_mode", strdatatype, att_space).write(strdatatype,
																		   strwritebuf);
	}
}


template<int dim>
inline void evolution_history(H5::H5File &file,
						  const typename history::EventAndMacro<dim> &event_history)
{
	std::string path = "/"+event_history.name;
	if(not H5Lexists(file.getId(), path.c_str(), H5P_DEFAULT))
		file.createGroup(path.c_str());

	{   /* --------- write macro evolution ----------- */

		using DataRow = history::MacroSummaryRow;

		H5::CompType mtype(sizeof(DataRow));
		mtype.insertMember("time", HOFFSET(DataRow, time), H5::PredType::NATIVE_DOUBLE);
		mtype.insertMember("total_strain", HOFFSET(DataRow, total_strain), H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("ext_stress", HOFFSET(DataRow, ext_stress), H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("av_vm_plastic_strain", HOFFSET(DataRow, av_vm_plastic_strain),
						   H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("std_vm_plastic_strain", HOFFSET(DataRow, std_vm_plastic_strain),
						   H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("av_vm_stress", HOFFSET(DataRow, av_vm_stress), H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("std_vm_stress", HOFFSET(DataRow, std_vm_stress),
						   H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("av_pressure", HOFFSET(DataRow, av_pressure), H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("std_pressure", HOFFSET(DataRow, std_pressure), H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("av_potential_energy", HOFFSET(DataRow, av_potential_energy),
						   H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("std_potential_energy", HOFFSET(DataRow, std_potential_energy),
						   H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("av_slip_threshold", HOFFSET(DataRow, av_slip_threshold),
						   H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("std_slip_threshold", HOFFSET(DataRow, std_slip_threshold),
						   H5::PredType::NATIVE_FLOAT);

		hsize_t d[] = {event_history.macro_evolution.size()};
		H5::DataSpace space(1, d);
		H5::DataSet dataset = file.createDataSet(path + "/macro_evolution", mtype, space);

		dataset.write(event_history.macro_evolution.data(), mtype);
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

		hsize_t d[] = {event_history.driving.size()};
		H5::DataSpace space(1, d);
		H5::DataSet dataset = file.createDataSet(path + "/driving_events", mtype, space);

		dataset.write(event_history.driving.data(), mtype);
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
		mtype.insertMember("activation_protocol", HOFFSET(DataRow, activation_protocol),
						   H5::PredType::NATIVE_UINT);

		hsize_t d[] = {event_history.plastic.size()};
		H5::DataSpace space(1, d);
		H5::DataSet dataset = file.createDataSet(path + "/plastic_events", mtype, space);

		dataset.write(event_history.plastic.data(), mtype);
	}

	//	   {   /* --------- write renewal event history ----------- */
	//		  using DataRow = typename history::EventAndMacro<dim>::RenewSlipRow;
	//	      H5::CompType mtype( sizeof(DataRow) );
	//	      mtype.insertMember( "index", HOFFSET(DataRow, index), H5::PredType::NATIVE_UINT);
	//	      mtype.insertMember( "element", HOFFSET(DataRow, element), H5::PredType::NATIVE_UINT);
	//	      mtype.insertMember( "threshold", HOFFSET(DataRow, threshold), H5::PredType::NATIVE_FLOAT);
	//	      mtype.insertMember( "slip_angle", HOFFSET(DataRow, slip_angle), H5::PredType::NATIVE_FLOAT);
	//
	//	      hsize_t d[] = {event_history.renew.size()};
	//	      H5::DataSpace space( 1, d );
	//	      H5::DataSet dataset = file.createDataSet(path+"/renew", mtype, space);
	//
	//	      dataset.write( event_history.renew.data(), mtype );
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
		unsigned int n_probe = x.first;
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
		H5::DataSet dataset = file.createDataSet("/patch_info/" + std::to_string(n_probe), datatype,
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

	{   /* --------- write local probing snapshots ----------- */

		if(not H5Lexists(file.getId(), (path + "/local_probe").c_str(), H5P_DEFAULT))
			file.createGroup(path + "/local_probe");

		using DataRow = typename mepls::patches::PatchPropertiesSnapshot<dim>::DataRow;
		H5::CompType mtype(sizeof(DataRow));
		mtype.insertMember("ref_element", HOFFSET(DataRow, ref_element), H5::PredType::NATIVE_UINT);
		mtype.insertMember("theta", HOFFSET(DataRow, theta), H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("failed", HOFFSET(DataRow, failed), H5::PredType::NATIVE_UINT);
		mtype.insertMember("ss_00", HOFFSET(DataRow, ss_xx), H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("ss_11", HOFFSET(DataRow, ss_yy), H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("ss_01", HOFFSET(DataRow, ss_xy), H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("oi_00", HOFFSET(DataRow, oi_xx), H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("oi_11", HOFFSET(DataRow, oi_yy), H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("oi_01", HOFFSET(DataRow, oi_xy), H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("oi_eps", HOFFSET(DataRow, oi_eps), H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("ee_00", HOFFSET(DataRow, ee_xx), H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("ee_11", HOFFSET(DataRow, ee_yy), H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("ee_01", HOFFSET(DataRow, ee_xy), H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("ss_pe", HOFFSET(DataRow, ss_pe), H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("oi_pe", HOFFSET(DataRow, oi_pe), H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("ee_pe", HOFFSET(DataRow, ee_pe), H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("x", HOFFSET(DataRow, x), H5::PredType::NATIVE_FLOAT);
		mtype.insertMember("y", HOFFSET(DataRow, y), H5::PredType::NATIVE_FLOAT);


		unsigned int n = 0;
		for(auto &snapshot : patch_prop_snapshots)
		{
			hsize_t d[] = {snapshot.data.size()};
			H5::DataSpace space(1, d);
			H5::DataSet dataset = file
				.createDataSet(path + "/local_probe/" + std::to_string(n++), mtype, space);
			dataset.write(snapshot.data.data(), mtype);

			H5::DataSpace att_space(H5S_SCALAR);
			dataset.createAttribute("desired_target", H5::PredType::NATIVE_DOUBLE, att_space)
				   .write(H5::PredType::NATIVE_DOUBLE, &snapshot.desired_target);
			dataset.createAttribute("recorded_target", H5::PredType::NATIVE_DOUBLE, att_space)
				   .write(H5::PredType::NATIVE_DOUBLE, &snapshot.recorded_target);
			dataset.createAttribute("index", H5::PredType::NATIVE_UINT, att_space)
				   .write(H5::PredType::NATIVE_UINT, &snapshot.output_index);
			dataset.createAttribute("N_probe", H5::PredType::NATIVE_UINT, att_space)
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


inline std::string make_filename(const parameters::Standard &p)
{
	std::ostringstream file_descriptor_ostrg;
	std::ostringstream L;
	std::ostringstream coupling_constant;
	std::ostringstream elastic_properties;
	std::ostringstream k;
	std::ostringstream nslip;
	std::ostringstream lambda;
	std::ostringstream temperature;
	std::ostringstream seed;

	L << "L_" << p.sim.Nx << "x" << p.sim.Ny;
	coupling_constant << "-c_" << std::fixed << std::setprecision(6) << p.mat.coupling_constant;
	nslip << "-nslip_" << p.mat.n_slip_systems;
	elastic_properties << "-avG_" << std::fixed << std::setprecision(3) << p.mat.average_G;
	lambda << "-lambda_" << std::fixed << std::setprecision(4) << p.mat.lambda;
	k << "-k_" << std::fixed << std::setprecision(4) << p.mat.k;
	temperature << "-T" << std::fixed << std::setprecision(2) << p.mat.temperature;
	seed << "-seed_" << p.sim.seed;

	file_descriptor_ostrg << L.str() << lambda.str() << k.str() << nslip.str()
						  << elastic_properties.str() << coupling_constant.str()
						  << temperature.str() << seed.str();

	return file_descriptor_ostrg.str();
}

} // namespace write


namespace quench
{

template<int dim>
void make_eshelby_prestress(std::vector<element::Anisotropic<dim> *> &elements,
							const parameters::Standard &p,
							std::mt19937 &generator)
{
	std::vector<dealii::SymmetricTensor<dim, 2>> eigenstrain(p.sim.Nx * p.sim.Ny);

	std::normal_distribution<double> normal_dist_shear(0., p.mat.prestress_std_shear);
	std::normal_distribution<double> normal_dist_pressure(0., p.mat.prestress_std_pressure);

	dealii::SymmetricTensor<dim, 2> eigenstrain_shear_0;
	dealii::SymmetricTensor<dim, 2> eigenstrain_shear_1;
	dealii::SymmetricTensor<dim, 2> eigenstrain_pressure;
	dealii::SymmetricTensor<dim, 2> eigenstrain_;

	for(unsigned int n = 0; n < eigenstrain.size(); ++n)
	{
		eigenstrain_shear_0[0][1] = normal_dist_shear(generator);

		eigenstrain_shear_1[0][0] = normal_dist_shear(generator);
		eigenstrain_shear_1[1][1] = -eigenstrain_shear_1[0][0];

		eigenstrain_pressure[0][0] = normal_dist_pressure(generator);
		eigenstrain_pressure[1][1] = eigenstrain_pressure[0][0];

		eigenstrain[n] = eigenstrain_shear_0 + eigenstrain_shear_1 + eigenstrain_pressure;
	}

	/* ---------------- add eigenstrain and get stress -----------------*/
	mepls::elasticity_solver::LeesEdwards<dim> solver(p.sim.Nx, p.sim.Ny);

	// since we only want to get a stress field with the desired values and that fulfills stress
	// equilibrium, it doesn't matter which elastic properties we used to obtain it. It will only
	// affect the values of the eigenstrain necessary to generate such stress field.
	assert(solver.get_n_elements() == elements.size());
	for(auto &element : elements)
		solver.set_elastic_properties(element->number(), element->C());

	solver.setup_and_assembly();

	for(unsigned int n = 0; n < eigenstrain.size(); ++n)
		solver.add_eigenstrain(n, eigenstrain[n]);

	solver.solve();


	// ensure that the average stress is 0 (it is not 0 because the eigenstrain is distributed with
	// average 0 but it fluctuates around it)
	dealii::SymmetricTensor<dim, 2> av_stress;
	auto stress = solver.get_stress();
	for(auto &tensor : stress)
		av_stress += tensor;
	av_stress /= double(stress.size());

	for(auto &tensor : stress)
		tensor -= av_stress;

	assert(stress.size() == elements.size());
	for(unsigned int n = 0; n < stress.size(); ++n)
		elements[n]->prestress(stress[n]);
}


template<int dim>
void equilibrate_initial_structure_rejection(mepls::system::System<dim> &system,
											 const parameters::Standard &p,
											 std::mt19937 &generator)
{
	// this ensures that initially all the thresholds are above the pre-stresses

	if(p.out.verbosity and omp_get_thread_num() == 0)
		std::cout << ">> Equilibrating initial structure... " << std::endl;

	mepls::element::RenewInstruct<dim> renew_instuct;

	for(auto &element : system)
	{
		bool unstable = true;
		while(unstable)
		{
			element->renew_structural_properties();

			unstable = false;
			for(auto &slip : *element)
				if(slip->barrier <= 0)
				{
					unstable = true;
					break;
				}
		}
	}

	for(auto &element : system)
		for(auto &slip : *element)
			assert(slip->barrier > 0);


	if(p.out.verbosity and omp_get_thread_num() == 0)
		std::cout << ">>>> Done " << std::endl;
}


template<int dim>
void equilibrate_initial_structure_relaxation(mepls::system::System<dim> &system,
											  const parameters::Standard &p)
{
	// this ensures that initially all the thresholds are above the pre-stresses

	if(p.out.verbosity and omp_get_thread_num() == 0)
		std::cout << ">> Equilibrating initial structure... " << std::endl;

	// relaxation dynamics to simulate an avalanche in which plastic deformation can occur.
	// This will modify the stress fields and a stable configuration might be eventually found.
	// Changes induced in the prestress should be accepatable since the prestress itself is
	// created from eigenstrains
	mepls::utils::ContinueSimulation continue_dummy;

	mepls::dynamics::relaxation(system, p.sim.fracture_limit, continue_dummy);M_Assert(
		continue_dummy(), "system reaches failure during equilibration");

	// as a consequente of the relaxation, now the elements have elastic strain fields giving rise
	// to stress. However, we don't want elastic fields when initiating the simulation. Therefore,
	// we convert the total local stress into a new prestress, and clean the rest of the
	// deformation history of the element
	for(auto &element : system)
	{
		dealii::SymmetricTensor<2, dim> new_prestress = element->stress();
		element->set_zero_deformation();
		element->prestress(new_prestress);
	}
	system.solver.clear();

	if(p.out.verbosity and omp_get_thread_num() == 0)
		std::cout << ">>>> Done " << std::endl;
}


template<int dim>
void run_thermal_evolution(mepls::system::System<dim> &system,
								history::EventAndMacro<dim> &history,
												const parameters::Standard &p,
												mepls::utils::ContinueSimulation &continue_simulation)
{
	mepls::dynamics::KMC<dim> kmc;
	auto &KMC_macro_evolution = history.macro_evolution;

	std::vector<double> rolling_av_stress;
	double rolling_av_stress_old = 0.;
	double rolling_av_stress_new = 0.;
	bool continue_kmc = true;
	unsigned int i = 0;
	unsigned int n = 100;

	while(continue_kmc)
	{
		++i;

		kmc(system, continue_simulation);
		mepls::dynamics::relaxation(system, p.sim.fracture_limit, continue_simulation);

		if(p.out.verbosity and omp_get_thread_num() == 0)
			std::cout << i << std::endl;

		history.add_macro(system);

		if(i % n == 0 and i > 1000)
		{
			rolling_av_stress_new = 0.;
			for(unsigned int j = KMC_macro_evolution.size() - n; j < KMC_macro_evolution.size();
				++j)
				rolling_av_stress_new += KMC_macro_evolution[j].av_vm_stress;
			rolling_av_stress_new /= double(n);

			continue_kmc = std::abs(
				rolling_av_stress_new - rolling_av_stress_old) / rolling_av_stress_old > 0.01;

			rolling_av_stress_old = rolling_av_stress_new;
		}
	}
}


template<int dim>
void convert_state_to_quench(mepls::system::System<dim> &system)
{
	// as a consequente of the evolution, now the elements have elastic strain fields giving rise
	// to stress. However, we don't want elastic fields when initiating the simulation. Therefore,
	// we convert the total local stress into a new prestress, and clean the rest of the
	// deformation history of the element and of the system
	for(auto &element : system)
	{
		dealii::SymmetricTensor<2, dim> new_prestress = element->stress();
		element->set_zero_deformation();
		element->prestress(new_prestress);
	}

	system.solver.clear();
	system.macrostate.clear();
}

} // namespace quench




template<int dim>
std::vector<element::Anisotropic<dim> *> create_elements(const parameters::Standard &p,
														 std::mt19937 &generator)
{
	auto weibull = [&](double x, double k, double lambda)
	{ return std::pow(x, k - 1) * std::exp(-std::pow(x / lambda, k)); };
	auto func = std::bind(weibull, std::placeholders::_1, p.mat.k_quench, p.mat.lambda_quench);
	auto threshold_distribution_ptr = mepls::utils::rand::create_distribution(func, 1e-3, 0., 1e-7,
																			  8);

	std::vector<element::Anisotropic<dim> *> elements;

	for(double n = 0; n < p.sim.Nx * p.sim.Ny; ++n)
	{
		typename element::Anisotropic<dim>::Config conf;

		conf.average_G = p.mat.average_G_quench;
		conf.average_K = p.mat.average_K;
		conf.weibull_shape_G = p.mat.weibull_shape_G;
		conf.weibull_shape_K = p.mat.weibull_shape_K;
		conf.alpha_tau = p.mat.alpha_tau;
		conf.n_slip_systems = p.mat.n_slip_systems_quench;
		conf.coupling_constant = p.mat.coupling_constant;
		conf.temperature = p.mat.temperature;
		conf.activation_rate = p.mat.activation_rate;
		conf.number = n;

		elements.push_back(
			new element::Anisotropic<dim>(*threshold_distribution_ptr, generator, conf));
	}

	if(p.mat.prestress)
		quench::make_eshelby_prestress<dim>(elements, p, generator);

	// renew the properties, so they take into account the prestress
	for(auto &element : elements)
		element->renew_structural_properties();


	return elements;
}


template<int dim>
void perform_reloading(mepls::system::System<dim> &system,
					   history::EventAndMacro<dim> &event_history,
					   bool is_forward,
					   const parameters::Standard &p)
{
	// Simulate using a copy of the original system; Clear the solver, and use the element stress
	// as prestress. Clean the deformation of the elements;
	system.solver.clear();

	mepls::element::Vector<dim> elements_replica;
	for(auto &element : system)
	{
		auto element_copy = element->make_copy();
		element_copy->set_zero_deformation();
		element_copy->prestress(element->stress());
		elements_replica.push_back(element_copy);
	}

	auto system_replica = system.get_new_instance(elements_replica, system.solver,
												  system.generator);
	auto &macrostate = system_replica->macrostate;
	system_replica->set_history(event_history);

	// initiate dynamics using copied system
	mepls::utils::ContinueSimulation continue_loading;
	while(continue_loading())
	{
		if(p.out.verbosity and omp_get_thread_num() == 0)
			std::cout << event_history.index << " | " << std::fixed << macrostate["total_strain"]
					  << " " << macrostate["ext_stress"] << " " << macrostate["pressure"]
					  << std::endl;

		event_history.add_macro( *system_replica );

		mepls::dynamics::finite_extremal_dynamics_step(1e-4 * 0.5, *system_replica, is_forward);
		mepls::dynamics::relaxation(*system_replica, p.sim.fracture_limit, continue_loading);

		continue_loading(std::abs(macrostate["total_strain"]) < 0.4 / 2., "System unloaded");
	}

	delete system_replica;

	if(p.out.verbosity and omp_get_thread_num() == 0)
	{
		std::cout << continue_loading << std::endl;
		std::cout << "Reloading " << (is_forward ? "forward" : "backward") << " finished"
				  << std::endl;
	}
}


} // namespace espci


#endif //_ESPCI_H

