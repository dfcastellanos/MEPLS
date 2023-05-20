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
#include <random>

int main()
{
    // the problem is in 2D
	constexpr unsigned int dim = 2;

    // the random number generator
    std::mt19937 generator(1234567);

	// the class of element we are using takes a configuration struct
	// For the moment, we will keep the default values
	example::element::Scalar<dim>::Config conf;

	// let's create the test element object
	example::element::Scalar<dim> element(conf, generator);

	//	A dealII's rank-2 symmsetric tensor. In this case, the stress tensor
	dealii::SymmetricTensor<2,dim> stress;
	stress[0][0] = 1.; // xx-component
	stress[1][1] = -1.; // yy-component
	stress[0][1] = 0.5; // xy-component

    // We set that stress in the element (a copy of the tensor will be stored)
	element.elastic_stress(stress);

    // We can ask in a similar way for its current stress tensor
	std::cout << "Stress = " << element.elastic_stress() << "\n\n";

	// The element's elastic properties are defined by a dealII's rank-4 stiffness tensor.
	// We can use a MEPLS function to easily create an isotropic one
	double G = 1.; // shear modulus
	double nu = 0.3; // poisson ratio
	dealii::SymmetricTensor<4,dim> C = mepls::utils::tensor::make_isotropic_stiffness<dim>(G, nu);

	element.C( C );

	// The element automatically updates its state according to linear elasticity.
	// Thus, we can ask e.g. for its elastic energy
	std::cout << "Elastic energy = " << element.energy_el() << "\n\n";

	//  We can iterate over the slips present in the element
	std::cout << "Slip angles = ";
	for(auto &slip : element)
		std::cout << slip->angle << ", ";
	std::cout << "\n\n";

	// We can ask for the number of slip objects
	std::cout << "Number of slip systems = " << element.size() << "\n\n";

	// and we can get a pointer to a specific slip as
	unsigned int slip_number = 1;
	mepls::slip::Slip<dim> * slip = element.slip( slip_number );

	std::cout << "Properties of slip system #" << slip_number << ": ";
    std::cout << "\n    * angle = " << slip->angle
			  << "\n    * threshold = " << slip->threshold;

	//  The slips have access to their parent's state, therefore the shear stress in the slip
	//  plane is already set
	std::cout << "\n    * shear stress = " << slip->eff_shear_stress;
  	std::cout << "\n    * barrier = " << slip->barrier << "\n\n";

	element.renew_structural_properties();

	// the previous slip doesn't exist anymore
	slip = element.slip( slip_number );

	std::cout << "Properties of the *new* slip system #" << slip_number << ": "
    		  << "\n    * angle = " << slip->angle
			  << "\n    * threshold = " << slip->threshold
			  << "\n    * shear stress = " << slip->eff_shear_stress
  			  << "\n    * barrier = " << slip->barrier << "\n\n";
}