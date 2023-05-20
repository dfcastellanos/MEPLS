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

#ifndef _EXAMPLE_H
#define _EXAMPLE_H


#include <mepls/element.h>

namespace example
{


namespace slip
{


template<int dim>
class Scalar : public mepls::slip::Slip<dim>
{
public:

	Scalar(double angle_, double threshold_)
		:
		mepls::slip::Slip<dim>()
	{
		threshold = threshold_;
		M_Assert(threshold > 0., "Expected threshold > 0");
		M_Assert(not std::isnan(threshold), "Threshold is NaN");

		angle = mepls::utils::mod(angle_, M_PI);
		M = mepls::utils::tensor::make_schmid<dim>(angle);
	};


	void update() override
	{
		eff_shear_stress = M * parent->stress();
		M_Assert(not std::isnan(eff_shear_stress), "eff_shear_stress is NaN");

		barrier = threshold - eff_shear_stress;
		activation_rate = std::exp(-barrier / temperature);
	}


	double get_critical_load_increment() override
	{
		M_Assert(barrier >= 0, "Barrier cannot be negative when computing the critical load increment.");

		return barrier / (M * parent->ext_stress_coeff());
	}

	dealii::SymmetricTensor<2, dim> get_eigenstrain_increment() override
	{
		return gamma * M;
	}

	using mepls::slip::Slip<dim>::angle;
	using mepls::slip::Slip<dim>::eff_shear_stress;
	using mepls::slip::Slip<dim>::pressure;
	using mepls::slip::Slip<dim>::threshold;
	using mepls::slip::Slip<dim>::barrier;
	using mepls::slip::Slip<dim>::parent;
	using mepls::slip::Slip<dim>::activation_rate;

	dealii::SymmetricTensor<2, dim> M;
	double gamma = 0.05;
	double temperature = 0.2;
};

} // namespace slip


namespace element
{

template<int dim>
class Scalar : public mepls::element::Element<dim>
{
public:

	struct Config
	{
		double gamma = 0.05;
		double temperature = 0.2;
		double k = 2;
		double lambda = 1.;
		unsigned int number = 0;
	};

	Scalar(Config &conf_, std::mt19937 &generator_)
		:
		mepls::element::Element<dim>(),
		generator(generator_),
		unif_distribution(0, 1),
		conf(conf_)
	{
		this->number(conf.number);

		mepls::element::Element<dim>::renew_structural_properties();
	}

	void renew_structural_properties_impl(mepls::event::Plastic<dim> &plastic_event) override
	{
		this->remove_slip_systems();

		double threshold = mepls::utils::rand::get_weibull_rand(conf.k, conf.lambda, unif_distribution(generator));

		auto slip1 = new slip::Scalar<dim>(0, threshold);
		slip1->temperature = conf.temperature;
		slip1->gamma = conf.gamma;
		this->add_slip_system(slip1);

		auto slip2 = new slip::Scalar<dim>(M_PI / 2., threshold);
		slip2->temperature = conf.temperature;
		slip2->gamma = conf.gamma;
		this->add_slip_system(slip2);
	}

	Config conf;

  private:
	std::mt19937 &generator;
	std::uniform_real_distribution<double> unif_distribution;
};

} // namespace element


} // namespace example


#endif //_EXAMPLE_H

