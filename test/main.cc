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

#include <mepls/utils.h>
#include <mepls/solver.h>

#include <iostream>

namespace mepls
{

namespace test
{

template<int dim>
double hydrostatic_eigenstrain_test(mepls::elasticity_solver::Solver<dim> &solver)
{
	solver.clear();

	dealii::SymmetricTensor<2,dim> eigenstrain;
	eigenstrain[0][0] = 1.;
	eigenstrain[1][1] = 1.;

	for(unsigned int n=0; n < solver.get_n_elements(); ++n)
		solver.add_eigenstrain(n, eigenstrain);

	solver.solve();

	auto &stress = solver.get_stress();
	dealii::SymmetricTensor<dim,2> S;
	for(unsigned int n=0; n < stress.size(); ++n)
		S += stress[n];
	S /= double(stress.size());

	double preassure = -(S[0][0]+S[1][1])/2.;

	return preassure;
}

template<int dim>
double dev_eigenstrain_test(mepls::elasticity_solver::Solver<dim> &solver)
{
	solver.clear();

	dealii::SymmetricTensor<2,dim> eigenstrain;
	eigenstrain[0][1] = 1.;

	for(unsigned int n=0; n < solver.get_n_elements(); ++n)
		solver.add_eigenstrain(n, eigenstrain);

	solver.solve();

	auto &stress = solver.get_stress();
	dealii::SymmetricTensor<dim,2> S;
	for(unsigned int n=0; n < stress.size(); ++n)
		S += stress[n];
	S /= double(stress.size());

	return S[0][1];
}

void shear_boundary_solver(bool write_vtu=false)
{
	std::cout << "\n======= TESTING SHEAR BOUNDARY SOLVER =========" <<	std::endl;


	/* ------ crete solver ------ */
	constexpr unsigned int dim = 2;
	unsigned int Nx = 30;
	unsigned int Ny = 30;

	mepls::elasticity_solver::ShearBoundary<dim> solver(Nx,Ny);
	solver.theta = 0.;

	/* ------ set elastic prop. ------ */
	double G = 25;
	double nu = 0.3;
	std::mt19937 generator(1234345);
	std::uniform_real_distribution<double> unif_dist(0,1);

	for(unsigned int n=0; n < solver.get_n_elements(); ++n)
		solver.set_elastic_properties(n,  utils::tensor::make_isotropic_stiffness<dim>(G,nu));

	solver.setup_and_assembly();

	double shear_stress = dev_eigenstrain_test(solver);
	double expected_shear_stress = -G * 2 * 1;
	assert_result("Shear stress from shear eigenstrain", shear_stress, expected_shear_stress, 1e-10);

	double preassure = hydrostatic_eigenstrain_test(solver);
	double K = 2*G*(1+nu)/(1-nu); // using the isotropic C defined in namespace utils, and considering that K := tr[Stress]/tr[eigenstrain]
	double expected_preassure = K;
	assert_result("Preassure from hyd. eigenstrain", preassure, expected_preassure, 1e-10);

	solver.clear();
	double epsxy = 1.;
	solver.add_load_increment(epsxy);
	solver.solve();
	double ext_stress = solver.get_external_stress();
	double ext_strain = solver.get_total_strain();
	double expected_ext_stress = 2*G*epsxy;
	double expected_ext_strain = epsxy;
	assert_result("External (loading) stress", ext_stress, expected_ext_stress, 1e-10);
	assert_result("External (loading) strain", ext_strain, expected_ext_strain, 1e-10);


	if(write_vtu)
		solver.write_vtu("./shear_boundary_solver.vtu");
}


} // namespace test

} // namespace mepls


int main()
{
	std::cout << std::endl;
	std::cout << "============================================" << std::endl;
	std::cout << "================ MEPLS TESTS ===============" << std::endl;
	std::cout << "============================================" << std::endl;

	mepls::test::shear_boundary_solver();

	std::cout << "============================================" << std::endl;
	std::cout << "=============== TESTS FINISHED =============" << std::endl;
	std::cout << "============================================" << std::endl;

	return 0;
}