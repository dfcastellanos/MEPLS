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

#ifndef _example_H
#define _example_H

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

/*! This namespace contains classes and functions to perform simulations based on example models. */
namespace example
{

namespace parameters
{

struct Material
{
	// threshold distribution
	double threshold_scale = 1.;
	double threshold_scale_quench = 1.5;
	double threshold_shape = 3.;

	// modulus distribution
	double shear_modulus = 30.;
	double poissons_ratio = 0.3;
	
	// others
	double coupling_constant = 1.;
	double temperature = 0.2;
	double activation_rate = 1.;

	void declare_entries(dealii::ParameterHandler &prm)
	{
		/*! Declare the variables of this struct into the input parser object. */

		prm.enter_subsection("Material");


		prm.declare_entry("shear_modulus", mepls::utils::str::to_string(shear_modulus),
						  dealii::Patterns::Double(), "");
		prm.declare_entry("poissons_ratio", mepls::utils::str::to_string(poissons_ratio),
						  dealii::Patterns::Double(), "");

		prm.declare_entry("threshold_scale", mepls::utils::str::to_string(threshold_scale),
						  dealii::Patterns::Double(), "");
		prm.declare_entry("threshold_scale_quench", mepls::utils::str::to_string(threshold_scale_quench),
						  dealii::Patterns::Double(), "");

		prm.declare_entry("threshold_shape", mepls::utils::str::to_string(threshold_shape), dealii::Patterns::Double(), "");

		prm.declare_entry("coupling_constant", mepls::utils::str::to_string(coupling_constant),
						  dealii::Patterns::Double(), "");
		prm.declare_entry("temperature", mepls::utils::str::to_string(temperature),
						  dealii::Patterns::Double(0.0), "");
		prm.declare_entry("activation_rate", mepls::utils::str::to_string(activation_rate),
						  dealii::Patterns::Double(0.0), "");


		prm.leave_subsection();

	}

	void load(dealii::ParameterHandler &prm)
	{
		/*! Load values parsed by the input parser object into the structs variables. */

		prm.enter_subsection("Material");


		shear_modulus = prm.get_double("shear_modulus");
		poissons_ratio = prm.get_double("poissons_ratio");

		threshold_scale = prm.get_double("threshold_scale");
		threshold_scale_quench = prm.get_double("threshold_scale_quench");
		threshold_shape = prm.get_double("threshold_shape");

		coupling_constant = prm.get_double("coupling_constant");
		temperature = prm.get_double("temperature");
		activation_rate = prm.get_double("activation_rate");


		prm.leave_subsection();
	}

};


/*! Struct with parameters related to the simulation setup. */
struct Simulation
{

	unsigned int n_rep = 1;

	unsigned int Nx = 32;
	/*!<  Number of elements in the x-direction of the rectangular lattice forming the system. */

	unsigned int Ny = 32;
	/*!< Number of elements in the y-direction of the rectangular lattice forming the system.  */

	unsigned int seed = 1234567;
	/*!<  Seed to initialized the random number engine. */

	double monitor_limit = 0.1;
	/*!<  Value of the magnitude defined by \ref monitor_name at which the simulation must stop. */

	double fracture_limit = 2;
	/*!<  Number of plastic events triggered during a call to \ref mepls::dynamics::relaxation
	 * divided by the number of elements composing the system at which the simulation must stop.
	 * A default value of 2 means that during a relaxation, every element has deformed twice,
	 * which is very big. Such value can, in most cases, indicate that the material has fractured
	 *  and therefore the simulation must stop. */


	void declare_entries(dealii::ParameterHandler &prm)
	{
		/*! Declare the variables of this struct into the input parser object. */
		prm.enter_subsection("Simulation Setup");

		prm.declare_entry("n_rep", mepls::utils::str::to_string(n_rep),
						  dealii::Patterns::Integer(0), "");
		prm.declare_entry("Nx", mepls::utils::str::to_string(Nx), dealii::Patterns::Integer(0), "");
		prm.declare_entry("Ny", mepls::utils::str::to_string(Ny), dealii::Patterns::Integer(0), "");
		prm.declare_entry("seed", mepls::utils::str::to_string(seed), dealii::Patterns::Integer(0),
						  "");
		prm.declare_entry("monitor_limit", mepls::utils::str::to_string(monitor_limit),
						  dealii::Patterns::Double(0.0), "");
		prm.declare_entry("fracture_limit", mepls::utils::str::to_string(fracture_limit),
						  dealii::Patterns::Double(0.0), "");


		prm.leave_subsection();
	}

	void load(dealii::ParameterHandler &prm)
	{
		/*! Load values parsed by the input parser object into the structs variables. */
		prm.enter_subsection("Simulation Setup");

		n_rep = prm.get_integer("n_rep");
		Nx = prm.get_integer("Nx");
		Ny = prm.get_integer("Ny");
		seed = prm.get_integer("seed");
		monitor_limit = prm.get_double("monitor_limit");
		fracture_limit = prm.get_double("fracture_limit");

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
class Example: public mepls::slip::Slip<dim>
{
public:

	Example(double angle_, double coupling_constant_, double temperature_)
		:
		mepls::slip::Slip<dim>()
	{
		angle = mepls::utils::mod(angle_, M_PI);
		coupling_constant = coupling_constant_;
		temperature = temperature_;

		M = mepls::utils::tensor::make_schmid<dim>(angle);
	}

	void update() override
	{
		eff_shear_stress = M * parent->stress();
		barrier = threshold - eff_shear_stress;
		activation_rate = activation_rate * std::exp(-barrier / temperature);

		M_Assert(threshold > 0., "Expected threshold > 0");
		M_Assert(not std::isnan(threshold), "Threshold is NaN");
		M_Assert(not std::isnan(eff_shear_stress), "eff_shear_stress is NaN");
	}

	double get_critical_load_increment() override
	{
		M_Assert(barrier >= 0,
				 "Negative barrier found while looking for critical load increment.");M_Assert(
			(parent->ext_stress_coeff()[0][0] != 0.) or (parent->ext_stress_coeff()[1][1] != 0.) or (parent->ext_stress_coeff()[0][1] != 0.),
			"Ext. stress coeffs. are all 0");

		return barrier / (M * parent->ext_stress_coeff() );
	}

	dealii::SymmetricTensor<2, dim> get_eigenstrain_increment() override
	{
		double eff_shear_stress_variation = - coupling_constant * eff_shear_stress;

		double A = M * parent->S() * M;
		double gamma = eff_shear_stress_variation / A; // limit drop is local shear stress
		dealii::SymmetricTensor<2, dim> eigenstrain_dev = gamma * M;

		return eigenstrain_dev;
	}

	using mepls::slip::Slip<dim>::angle;
	using mepls::slip::Slip<dim>::eff_shear_stress;
	using mepls::slip::Slip<dim>::threshold;
	using mepls::slip::Slip<dim>::barrier;
	using mepls::slip::Slip<dim>::parent;
	using mepls::slip::Slip<dim>::activation_rate;

	dealii::SymmetricTensor<2, dim> M;
	double coupling_constant = 0.;
	double temperature = 0.2;
};

} // namespace slip


namespace element
{

template<int dim>
class Example: public mepls::element::Element<dim>
{
public:

	struct Config
	{
		double temperature = 1.;
		double activation_rate = 1.;
		double coupling_constant = 1.;
		double threshold_shape = 1.;
		double threshold_scale = 1.;
		double threshold_scale_quench = 1.;
		unsigned int number = 0;
		dealii::SymmetricTensor<4,dim> C;
	};


	Example(std::mt19937 &generator_,
				const Config &conf_)
		:
		mepls::element::Element<dim>(),
		generator(generator_),
		unif_distribution(0, 1),
		conf(conf_)
	{
		this->number(conf.number);
		this->C(conf.C);

		auto slip_system_1 = new slip::Example<dim>(0, conf.coupling_constant, conf.temperature);
		slip_system_1->threshold = mepls::utils::rand::get_weibull_rand(conf.threshold_shape, conf.threshold_scale_quench, unif_distribution(generator));

		auto slip_system_2 = new slip::Example<dim>(M_PI / 2., conf.coupling_constant, conf.temperature);
		slip_system_2->threshold = mepls::utils::rand::get_weibull_rand(conf.threshold_shape, conf.threshold_scale_quench, unif_distribution(generator));

		this->add_slip_system(slip_system_1);
		this->add_slip_system(slip_system_2);
	}


	void renew_structural_properties_impl(mepls::element::RenewInstruct<dim> &renew_instruct) override
	{
		for(auto &slip : *this)
		{
			slip->threshold = mepls::utils::rand::get_weibull_rand(conf.threshold_shape, conf.threshold_scale, unif_distribution(generator));
			slip->update();
		}
	}


  private:
	std::mt19937 &generator;
	std::uniform_real_distribution<double> unif_distribution;
	Config conf;
};

} // namespace element



namespace write
{

inline void file_attrs(H5::H5File &file, const parameters::Parameters &p)
{
	H5::DataSpace att_space(H5S_SCALAR);
	file.createAttribute("threshold_scale", H5::PredType::NATIVE_DOUBLE, att_space).write(
		H5::PredType::NATIVE_DOUBLE, &p.mat.threshold_scale);
	file.createAttribute("threshold_shape", H5::PredType::NATIVE_DOUBLE, att_space).write(
		H5::PredType::NATIVE_DOUBLE, &p.mat.threshold_shape);

	file.createAttribute("shear_modulus", H5::PredType::NATIVE_DOUBLE, att_space).write(
		H5::PredType::NATIVE_DOUBLE, &p.mat.shear_modulus);

	file.createAttribute("poissons_ratio", H5::PredType::NATIVE_DOUBLE, att_space).write(
		H5::PredType::NATIVE_DOUBLE, &p.mat.poissons_ratio);

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

}


template<int dim>
inline void evolution_history(H5::H5File &file,
						  const typename mepls::history::History<dim> &history)
{
	std::string path = "/"+history.name();
	if(not H5Lexists(file.getId(), path.c_str(), H5P_DEFAULT))
		file.createGroup(path.c_str());

	{   /* --------- write macro evolution ----------- */

		using DataRow = mepls::history::MacroSummaryRow;

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
		mtype.insertMember("av_potential_energy", HOFFSET(DataRow, av_potential_energy),
						   H5::PredType::NATIVE_DOUBLE);
		mtype.insertMember("std_potential_energy", HOFFSET(DataRow, std_potential_energy),
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

		hsize_t d[] = {history.macro_evolution.size()};
		H5::DataSpace space(1, d);
		H5::DataSet dataset = file.createDataSet(path + "/macro_evolution", mtype, space);

		dataset.write(history.macro_evolution.data(), mtype);
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
		mtype.insertMember("activation_protocol", HOFFSET(DataRow, activation_protocol),
						   H5::PredType::NATIVE_UINT);

		hsize_t d[] = {history.plastic.size()};
		H5::DataSpace space(1, d);
		H5::DataSet dataset = file.createDataSet(path + "/plastic_events", mtype, space);

		dataset.write(history.plastic.data(), mtype);
	}

}



template<int dim>
inline void snapshots(H5::H5File &file,
					  std::string path,
					  const std::vector<mepls::snapshot::Stress<dim> > &stress_snapshots)
{

	if(not H5Lexists(file.getId(), path.c_str(), H5P_DEFAULT))
		file.createGroup(path);


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


inline std::string make_filename(const parameters::Parameters &p)
{
	std::ostringstream file_descriptor_ostrg;
	std::ostringstream L;
	std::ostringstream coupling_constant;
	std::ostringstream elastic_properties;
	std::ostringstream threshold_shape;
	std::ostringstream threshold_scale;
	std::ostringstream temperature;
	std::ostringstream seed;

	L << "L_" << p.sim.Nx << "x" << p.sim.Ny;
	coupling_constant << "-c_" << std::fixed << std::setprecision(3) << p.mat.coupling_constant;
	elastic_properties << "-G_" << std::fixed << std::setprecision(3) << p.mat.shear_modulus;
	threshold_scale << "-k_" << std::fixed << std::setprecision(3) << p.mat.threshold_scale;
	threshold_shape << "-k_" << std::fixed << std::setprecision(3) << p.mat.threshold_shape;
	temperature << "-T_" << std::fixed << std::setprecision(3) << p.mat.temperature;
	seed << "-seed_" << p.sim.seed;

	file_descriptor_ostrg << L.str() << threshold_scale.str() << threshold_shape.str()
						  << elastic_properties.str() << coupling_constant.str()
						  << temperature.str() << seed.str();

	return file_descriptor_ostrg.str();
}

} // namespace write



} // namespace example


#endif //_example_H

