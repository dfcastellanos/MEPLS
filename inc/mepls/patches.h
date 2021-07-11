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

#ifndef __PATCHES_H_
#define __PATCHES_H_


#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/base/point.h>
#include <deal.II/base/timer.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>
#include <string>
#include <map>

#include <mepls/utils.h>
#include <mepls/solver_impl.h>
#include <mepls/dynamics.h>


namespace mepls
{


/*! This namespace contains tools for creating and testing patches. Patches are
 * connected sub-domains of the full system, which are cut-out and mechanically
 * shear-tested in isolation (see @cite DFCastellanos_CRP @cite Patinet2016 @cite Barbot2018).
 *  In elasto-plastic models, patches are composed by mesoscale elements and
 *  represented as a @ref system::System object, with the same properties and evolution laws as
 *  the full-size system from where they are cut-out. Performing shear tests on individual patches
 * provides statistical information about the spatially-fluctuating local
 * mechanical response of the system. While full systems can be periodic or
 * non-periodic patches are, by definition, non-periodic. Therefore, computing
 * elastic fields within patches require solving elastic problems in the
 * presence of surfaces.
 *
 * Three different states during a shear test are distinguished:
 * 1. Stable sate ("ss"). The patch state previous to the initiation of the 
 * test.
 * 2. The onset of instability ("oi"). The state of the patch at the very
 * moment of the first mechanical instability (i.e., at the moment of
 * activation of the first slip event within the patch).
 * 3. End-of-event state ("ee"). The state of the patch when the plastic
 * deformation initiated in the oi-state has finished. Note that such
 * plastic deformation can, in general, be composed of several slip
 * systems beyond the initial one triggered in the oi-state. */
namespace patches
{


/*! A struct containing the properties that are measured when a patch is 
 * tested. */
template<int dim>
struct PatchPropertiesTensorial
{
	unsigned int ref_element = 0;
	/*!<  Number of the element used as a reference when creating a patch. It can
	 * be used a patch id to distinguish it from other patches. */

	double theta = 0.;
	/*!< Orientation of the shear test performed on the patch. This angle is
	 * the same used in @ref utils::tensor::schmid. */

	dealii::SymmetricTensor<2, dim> stress_ss;
	/*!< Stress tensor averaged over the elements composing the patch. The
	 * operation is performed in the ss-state. */

	dealii::SymmetricTensor<2, dim> stress_oi;
	/*!< Stress tensor averaged over the elements composing the patch. The
	 * operation is performed in the oi-state. */

	dealii::SymmetricTensor<2, dim> stress_ee;
	/*!< Stress tensor averaged over the elements composing the patch. The
	 * operation is performed in the ee-state. */

	double resolved_elastic_shear_strain_oi;
	/*!< Elastic shear strain applied to the patch, with shear orientation
	 * @ref theta, to reach the oi-state from the initial ss-state. */

	double energy_el_ss = 0.;
	/*!< Elastic energy averaged over the elements composing the patch. The
	 * operation is performed in the ss-state. */

	double energy_el_oi = 0.;
	/*!< Elastic energy averaged over the elements composing the patch. The
	 * operation is performed in the oi-state. */

	double energy_el_ee = 0.;
	/*!< Elastic energy averaged over the elements composing the patch. The
	 * operation is performed in the ee-state. */

	double energy_conf_ss = 0.;
	/*!< Configurational energy averaged over the elements composing the patch. The
	 * operation is performed in the ss-state. */

	double energy_conf_oi = 0.;
	/*!< Configurational energy averaged over the elements composing the patch. The
	 * operation is performed in the oi-state. */

	double energy_conf_ee = 0.;
	/*!< Configurational energy averaged over the elements composing the patch. The
	 * operation is performed in the ee-state. */

	std::vector<double> coords;
	/*!< Cartesian coordinates of the center of the patch. */

	bool failed = false;
	/*!< Whether the patch has undergone mechanical failure during the shear
	 * test. This can be the case if, e.g., strain-softening is present in the
	 * system. */

};


template<int dim>
inline std::vector<std::pair<int, int> > get_relative_coordinates(int n)
{
	/*! Get the integer coordinates relative to the origin \f$ (0,0) \f$ up to
	 * an integer distance \f$ n \f$, given by the cartesian product
	 * \f$ [0,n] \times [0,n] \f$.
	 *
	 * <hr>
	 *
	 *  For example, for \f$ n=2 \f$,
	 *
	 * @code
	  ----------------------------
	  | (0, 2) | (1, 2) | (2, 2) |
	  ----------------------------
	  | (0, 1) | (1, 1) | (2, 1) |
	  ----------------------------
	  | (0, 0) | (1, 0) | (2, 0) |
	  ----------------------------
	 * @endcode
	 *
	 */

	assert(dim == 2);M_Assert(n >= 1, "Patch size must be bigger or equal to 1");

	std::vector<std::pair<int, int> > rel_neighbor_coords;

	for(int i = 0; i < n; ++i)
		for(int j = 0; j < n; ++j)
			rel_neighbor_coords.push_back(std::make_pair(i, j));

	return rel_neighbor_coords;
}


template<int dim>
inline void select_patch_elements(
	element::Vector<dim> &patch_elements,
	std::vector<double> &patch_center_coords,
	const element::Element<dim> *ref_element,
	const element::Vector<dim> &elements,
	unsigned int Nx,
	unsigned int Ny,
	unsigned int N)
{
	/*! Create a square patch by selecting the elements composing it.
	 *
	 * @param patch_elements elements composing the patch
	 * @param patch_center_coords cartesian coordinates of the patch center
	 * @param ref_element reference element
	 * @param elements elements composing the full system
	 * @param Nx number of elements in the full system along direction x
	 * @param Ny number of elements in the full system along direction y
	 * @param N number of elements per side (i.e. the total number of elements
	 * in the patch is \f$ N^{\rm dim} \f$)
	 *
	 * @warning the patches are created assuming that the full system is
	 * periodic.
	 *
	 * <hr>
	 *
	 * Example of the creation of a patch with N=3 and reference element 18 from
	 * a square 5x5 bi-periodic system with 25 elements (the numbers between
	 * stars denote the elements chosen for composing the patch):
	 *
	 * @code
	  --------------------------
	  |*20*| 21 | 22 |*23*|*24*|
	  --------------------------
	  |*15*| 16 | 17 |*18*|*19*|
	  --------------------------
	  | 10 | 11 | 12 | 13 | 14 |
	  --------------------------
	  |  5 |  6 |  7 |  8 |  9 |
	  --------------------------
	  | *0*|  1 |  2 | *3*| *4*|
	  --------------------------
	 * @endcode
	 *
	 * The coordinates of the elements relative to the reference element 18 are
	 * 18->[0,0], 24->[1,1] etc. The cut-out patch is then:
	 *
	 * @code
	  ----------------
	  |  3 |  4 |  0 |
	  ----------------
	  | 23 | 24 | 20 |
	  ----------------
	  | 18 | 19 | 15 |
	  ----------------
	 * @endcode
	 *
	 * The selected patch becomes a system in its own, with an element numbering
	 * following the same convention as the full system, that is:
	 *
	 * @code
	  --------------
	  | 6 |  7 | 8 |
	  --------------
	  | 3 |  4 | 5 |
	  --------------
	  | 0 |  1 | 2 |
	  --------------
	 * @endcode
	 *
	 * Thus, the vector of elements @ref system::System<dim>.elements in the
	 * patch system has the elements ordered as [18,19,15,23,24,20,3,4,0]
	 */

	std::vector<std::pair<int, int> > rel_neighbor_coords = get_relative_coordinates<dim>(N);

	patch_elements.clear();
	patch_elements.resize(rel_neighbor_coords.size());

	// coordinates of the reference element; map 1D index into 2D array
	int c_x = (ref_element->number()) % Nx;
	int c_y = std::floor((ref_element->number()) / Nx);

	// number of elements in the patch along x and y directions
	unsigned int nx = N;
	unsigned int ny = N;

	M_Assert(rel_neighbor_coords.size() == nx * ny, "Unexpected number of patch elements");

	// patch center coordinates
	double xm = 0.;
	double ym = 0.;

	for(auto &d : rel_neighbor_coords)
	{
		// coordinates of the neighbor relatives to the ref. element
		int dx = d.first;
		int dy = d.second;

		// coordinates of the neighbor in the full system
		int x_ = dx + c_x;
		int y_ = dy + c_y;

		// wrap the coordinates around the full system boundaries to consider
		// a periodic system
		int x = utils::mod(x_, Nx);
		int y = utils::mod(y_, Ny);

		// to compute the average of the (non-warped) coords. of the patch
		// elements. The average patch coords. is the geometrical center of the
		// patch
		xm += x_;
		ym += y_;

		// element number of the neighbor in the full system numbering (this
		// just maps 2D coordinates to a 1D array).
		unsigned int n_neighbor = x + Nx * y;

		M_Assert(elements[n_neighbor]->number() == n_neighbor,
				 "elements are not correctly ordered");

		// we need the elements in the orderd following the same convention as
		// in the full system. Following the example above, the vector
		// must have the elements in the order [18,19,15,23,24,20,3,4,0]

		// element number of the neighbor in the small patch-system numbering
		unsigned int n_neighbor_probing = dx + nx * dy;
		patch_elements[n_neighbor_probing] = elements[n_neighbor];
	}

	// wrap the average patch coordinates to account for the periodicity of the
	// full system
	patch_center_coords.resize(2);
	patch_center_coords[0] = utils::mod(xm / double(rel_neighbor_coords.size()), Nx);
	patch_center_coords[1] = utils::mod(ym / double(rel_neighbor_coords.size()), Ny);
}


template<int dim>
std::map<unsigned int, std::vector<std::vector<unsigned int>>> make_patch_to_element_map(
	const element::Vector<dim> &elements,
	const std::vector<int> &N_list,
	unsigned int Nx,
	unsigned int Ny)
{
	/*! Returns a map between patch size and all the patches of that size. The
	 * keys of the map correspond to the available patch sizes (which are
	 * defined in input vector N_list). The values of the map are vectors
	 * of vectors v[i][j], where i corresponds to the number of the reference
	 * element used for creating the patch and j to the number of the element.
	 *
	 * @param elements composing the full system
	 * @param N_list vector of patch sizes to include in the map
	 * @param Nx number of elements in the full system along direction x
	 * @param Ny number of elements in the full system along direction y
	 * */

	std::map<unsigned int, std::vector<std::vector<unsigned int>>> patch_to_element_map;
	for(auto N : N_list)
	{
		patch_to_element_map[N] = std::vector<std::vector<unsigned int>>();
		std::vector<std::pair<int, int> > patch_elements_rel_coords = get_relative_coordinates<dim>(
			N);

		for(auto &reference_element : elements)
		{
			element::Vector<dim> patch_elements;
			std::vector<double> patch_center_coords;
			select_patch_elements<dim>(patch_elements, patch_center_coords, reference_element,
									   elements, Nx, Ny, N);

			std::vector<unsigned int> patch_elements_numbers;
			for(auto &patch_element : patch_elements)
				patch_elements_numbers.push_back(patch_element->number());

			patch_to_element_map[N].push_back(patch_elements_numbers);

			// make sure that the index of the first vector corresponds to the
			// patch ref. element number
			assert(patch_to_element_map[N].size() == reference_element->number() + 1);
		}
	}

	return patch_to_element_map;
}


template<int dim>
inline void get_averages(const element::Vector<dim> &elements,
						dealii::SymmetricTensor<2, dim> &av_stress,
						double &av_energy_el,
						double &av_energy_conf)
{
	/*! Compute the average properties over the input elements. */

	av_stress.clear();
	av_energy_el = 0.;
	av_energy_conf = 0.;

	for(auto &element : elements)
	{
		av_stress += element->stress();
		av_energy_el += element->energy_el();
		av_energy_conf += element->energy_conf();
	}

	av_stress /= double(elements.size());
	av_energy_el /= double(elements.size());
	av_energy_conf /= double(elements.size());
}


template<int dim>
void apply_patch_shear_test(
	PatchPropertiesTensorial<dim> &patch_properties,
	system::System<dim> &patch_system,
	utils::ContinueSimulation &continue_shear_test,
	bool do_ee,
	double dgamma)
{
	/*! Apply a shear test to the input system. The system is sheared in discrete
	 * shear strain increments along with a certain orientation. The properties
	 * of the system at the ss- oi- and ee-states are computed.
	 *
	 * @param results of the shear test
	 * @param system on which we apply the shear test
	 * @param do_ee if false, compute only the ss- and oi-state. If true, compute
	 * also the ee-state
	 * @param dgamma amplitude of the discrete external shear strain increments
	 * @param Nx number of elements in the input system along direction x
	 * @param Ny number of elements in the input system along direction y
	 *
	 * @note this function is normally used on systems representing patches,
	 * but can be used, in general, on any system.
	 */

	auto &elements = patch_system.elements;
	auto &solver = patch_system.solver;

	// stable state
	get_averages(elements, patch_properties.stress_ss, patch_properties.energy_el_ss,
	             patch_properties.energy_conf_ss);

	// onset instability state: drive the system until it becomes unstable
	// (i.e., a slip system has reached its threshold).
	double deps = 0.5 * dgamma;
	dynamics::finite_extremal_dynamics_step(deps, patch_system);
	get_averages(elements, patch_properties.stress_oi, patch_properties.energy_el_oi,
	             patch_properties.energy_conf_oi);
	patch_properties.resolved_elastic_shear_strain_oi = solver.get_total_strain();

	if(do_ee)
	{
		// ee-state: we let the system relax by performing the slip event
		// triggered in the oi-state (and possibly, other events triggered
		// induced by it)
		dynamics::relaxation(patch_system, continue_shear_test);
	    get_averages(elements, patch_properties.stress_ee, patch_properties.energy_el_ee,
	                 patch_properties.energy_conf_ee);
	}
	else
	{
		// use oi as ee, so the variation from one to the other is zero
		patch_properties.stress_ee = patch_properties.stress_oi;
		patch_properties.energy_el_ee = patch_properties.energy_el_oi;
		patch_properties.energy_conf_ee = patch_properties.energy_conf_oi;
	}
}


template<int dim>
void analyze_patch_ensemble(
	std::vector<PatchPropertiesTensorial<dim>> &data_patch,
	const system::System<dim> &system,
	unsigned int N,
	const std::vector<double> &theta_list,
	bool do_ee)
{
	/*! Apply shear tests to an ensemble of patches of size N. A patch of
	 * size N is created for each reference element (see @ref
	 * select_patch_elements for the interpretation) considered. Reference
	 * elements are defined by iterating over all the elements in the system and
	 * selecting a reference element for every N/2 elements. In this way, we
	 * reduce the spatial overlap between different patches.
	 *
	 * @param data_patch vector containing the properties measured for each patch
	 * @param system from where the patches are created
	 * @param N patch size (see @ref select_patch_elements)
	 * @param theta_list vector with the shear orientations along which shear
	 * tests are to be performed
	 * @param do_ee if false, compute only the ss- and oi-state. If true, compute
	 * also the ee-state
	 */

	const auto &full_system_solver = system.solver;
	const unsigned int Nx = full_system_solver.get_Nx();
	const unsigned int Ny = full_system_solver.get_Ny();
	auto &generator = system.generator;
	auto &elements = system.elements;

	// pointers to elements composing the patch which is to be tested
	element::Vector<dim> patch_elements;

	// pointers to copies of the probing elements to avoid modifying the
	// originals
	element::Vector<dim> patch_elements_copy;
	std::vector<double> patch_center_coords;

	for(auto &theta : theta_list)
	{
		elasticity_solver::ShearBoundary<dim> probring_solver(N, N);
		probring_solver.theta = theta;
		probring_solver.setup_and_assembly();

		/* ------- make patches and test them ------- */

		// we select random non-repeated reference elements. We build a patch using the
		// reference element's neighborhood using select_patch_elements(). The number of
		// elements per patch is N^2, therefore number of elements in the system divided by N^2
		// gives a number of patches such that their overlap is, on average, not too big
		unsigned int n_patches = int( double(elements.size())/double(N*N) );
		assert(n_patches >= 1);

		mepls::element::Vector<dim> reference_elements;

   		 std::sample(elements.begin(), elements.end(), std::back_inserter(reference_elements),
   		             n_patches, generator);

		for(auto &reference_element : reference_elements)
		{
			select_patch_elements<dim>(patch_elements, patch_center_coords, reference_element,
									   elements, Nx, Ny, N);

			M_Assert(patch_elements.size() == probring_solver.get_n_elements(), "");

			/* ------- make copy of elements to work on them safely -------*/
			std::vector<unsigned int> original_element_numbers(patch_elements.size());
			patch_elements_copy.resize(patch_elements.size());

			for(unsigned int n = 0; n < patch_elements.size(); ++n)
			{
				auto element = patch_elements[n];
				auto element_copy = element->make_copy();

				element_copy->state_to_prestress();

				original_element_numbers[n] = element->number();

				// number it according to the patch system 1D index, which
				// correspond with the order in the vector patch_elements
				element_copy->number(n);
				patch_elements_copy[n] = element_copy;
			}

			// the patch system is different from the full system. We need a
			// specific solver for the patch system. However, we can reuse the
			// assembly of the full system since, except for the system
			// boundaries, the properties of both systems are locally the same
			probring_solver.copy_assembly(original_element_numbers, full_system_solver);

			element::calculate_ext_stress_coefficients(patch_elements_copy, probring_solver);

			// @ref element::Element<dim>.S_ needs to be computed only if plastic
			// deformation takes place. This is only the case if we change from
			// oi- to the ee-state
			if(do_ee)
				element::calculate_local_stress_coefficients(patch_elements_copy, probring_solver);


			/* -------- apply shear test to the patch --------*/
			system::System<dim> *probing_system = system.get_new_instance(patch_elements_copy,
																		  probring_solver,
																		  generator);
			utils::ContinueSimulation continue_shear_test;
			PatchPropertiesTensorial<dim> patch_properties;
			patch_properties.ref_element = reference_element->number();
			patch_properties.theta = theta * 180. / M_PI;
			patch_properties.coords = patch_center_coords;

			double dgamma = 1e-4;
			apply_patch_shear_test(patch_properties, *probing_system, continue_shear_test, do_ee,
								   dgamma);

#ifdef DEBUG
			dealii::SymmetricTensor<2,dim> M = utils::tensor::make_schmid<dim>(theta);
			double eff_ss = M*patch_properties.stress_ss;
			double eff_oi = M*patch_properties.stress_oi;
			M_Assert(eff_oi>eff_ss, "resolved shear stress on the shearing plane at the onset of instability is lower than at the stable state");
#endif

			if(not continue_shear_test())
			{
				patch_properties.failed = true;

#ifdef DEBUG
				std::cout << "* Warning: probing system failed" << std::endl;
#endif
			}

			data_patch.push_back(patch_properties);

			delete probing_system;

			for(auto &element : patch_elements_copy)
				delete element;

			patch_elements_copy.clear();

			probring_solver.clear();

		} // reference element loop
	} // theta loop

}


template<int dim>
void analyze_patch_ensemble_opt(
	std::vector<PatchPropertiesTensorial<dim>> &data_patch,
	const system::System<dim> &system,
	unsigned int N,
	const std::vector<double> &theta_list,
	bool do_ee)
{
	/*! It does the same as @ref analyze_patch_ensemble, but precalculates
	 * @ref element::Element<dim>.ext_stress_coeff_and @ref element::Element<dim>
	 * .S_ once and is reused for all the patches. This optimization is only
	 * valid under the same conditions discussed for @ref
	 * element::calculate_local_stress_coefficients_central, namely the full
	 * system has elastic homogeneous properties and periodic boundary
	 * conditions.
	 */


	// see @ref analyze_patch_ensemble for a more detailed documentation

	const auto &full_system_solver = system.solver;
	const unsigned int Nx = full_system_solver.get_Nx();
	const unsigned int Ny = full_system_solver.get_Ny();
	auto &generator = system.generator;
	auto &elements = system.elements;

	element::Vector<dim> patch_elements;
	element::Vector<dim> patch_elements_copy;
	std::vector<double> patch_center_coords;

	for(auto &theta : theta_list)
	{
		elasticity_solver::ShearBoundary<dim> probring_solver(N, N);
		probring_solver.theta = theta;
		probring_solver.setup_and_assembly();


		/* -------- precalculate once S and ext_stress_coeffs --------*/

		std::vector<dealii::SymmetricTensor<2, dim>> ext_stress_coefficients_vector;
		std::vector<dealii::SymmetricTensor<4, dim>> S_vector;
		{
			auto &reference_element = elements[0];
			select_patch_elements<dim>(patch_elements, patch_center_coords, reference_element,
									   elements, Nx, Ny, N);

			M_Assert(patch_elements.size() == probring_solver.get_n_elements(), "");

			std::vector<unsigned int> original_element_numbers(patch_elements.size());
			patch_elements_copy.resize(patch_elements.size());
			for(unsigned int n = 0; n < patch_elements.size(); ++n)
			{
				auto element = patch_elements[n];
				auto element_copy = element->make_copy();

				original_element_numbers[n] = element->number();

				element_copy->number(n);
				patch_elements_copy[n] = element_copy;
			}

			probring_solver.copy_assembly(original_element_numbers, full_system_solver);

			element::calculate_ext_stress_coefficients(patch_elements_copy, probring_solver);
			if(do_ee)
				element::calculate_local_stress_coefficients(patch_elements_copy, probring_solver);

			// extract the quantities of interest from them
			for(auto &element : patch_elements_copy)
			{
				ext_stress_coefficients_vector.push_back(element->ext_stress_coeff());
				S_vector.push_back(element->S());
				delete element;
			}

			patch_elements_copy.clear();

		}

		/* ------- proceed as in @ref analyze_patch_ensemble ------- */

		// we select random non-repeated reference elements. We build a patch using the
		// reference element's neighborhood using select_patch_elements(). The number of
		// elements per patch is N^2, therefore number of elements in the system divided by N^2
		// gives a number of patches such that their overlap is, on average, not too big
		unsigned int n_patches = int( double(elements.size())/double(N*N) );
		assert(n_patches >= 1);

		mepls::element::Vector<dim> reference_elements;

   		 std::sample(elements.begin(), elements.end(), std::back_inserter(reference_elements),
   		             n_patches, generator);

		for(auto &reference_element : reference_elements)
		{
			select_patch_elements<dim>(patch_elements, patch_center_coords, reference_element,
									   elements, Nx, Ny, N);

			M_Assert(patch_elements.size() == probring_solver.get_n_elements(), "");

			patch_elements_copy.resize(patch_elements.size());
			for(unsigned int n = 0; n < patch_elements.size(); ++n)
			{
				auto element = patch_elements[n];
				auto element_copy = element->make_copy();

				element_copy->state_to_prestress();

				element_copy->number(n);

				// use the values precalculated above

				element_copy->ext_stress_coeff(ext_stress_coefficients_vector[n]);
				element_copy->S(S_vector[n]);


				patch_elements_copy[n] = element_copy;
			}

			//		auto element_o = patch_elements_copy[0];
			//		auto slip_o = *(element_o->begin());
			//		std::cout << "ORIGINAL" << std::endl;
			//		std::cout << "prestress: " << element_o->prestress() << std::endl;
			//		std::cout << "elastic_strain: " << element_o->elastic_strain() << std::endl;
			//		std::cout << "stress: " << element_o->stress() << std::endl;
			//		std::cout << "slip: " << slip_o->angle << " " << slip_o->eff_shear_stress
			//				  << " " << slip_o->threshold << " " << slip_o->barrier
			//				  << std::endl << std::endl;
			//
			//		auto element = patch_elements_copy[0];
			//		auto slip = *(element->begin());
			//		std::cout << "BEFORE" << std::endl;
			//		std::cout << "prestress: " << element->prestress() << std::endl;
			//		std::cout << "elastic_strain: " << element->elastic_strain() << std::endl;
			//		std::cout << "stress: " << element->stress() << std::endl;
			//		std::cout << "slip: " << slip->angle << " " << slip->eff_shear_stress
			//				 << " " << slip->threshold << " " << slip->barrier << std::endl
			//				 << std::endl;

			/* -------- apply test --------*/
			system::System<dim> *probing_system = system.get_new_instance(patch_elements_copy,
																		  probring_solver,
																		  generator);
			utils::ContinueSimulation continue_shear_test;
			PatchPropertiesTensorial<dim> patch_properties;
			patch_properties.ref_element = reference_element->number();
			patch_properties.theta = theta * 180. / M_PI;
			patch_properties.coords = patch_center_coords;

			double dgamma = 1e-4;
			apply_patch_shear_test(patch_properties, *probing_system, continue_shear_test, do_ee,
								   dgamma);

			//		slip = *(element->begin());
			//		std::cout << "AFTER" << std::endl;
			//		std::cout << "prestress: " << element->prestress() << std::endl;
			//		std::cout << "elastic_strain: " << element->elastic_strain() << std::endl;
			//		std::cout << "stress: " << element->stress() << std::endl;
			//		std::cout << "slip: " << slip->angle << " " << slip->eff_shear_stress
			//				  << " " << slip->threshold << " " << slip->barrier << std::endl
			//				  << std::endl;
			//		abort();


#ifdef DEBUG
			dealii::SymmetricTensor<2,dim> M = utils::tensor::make_schmid<dim>(theta);
			double eff_ss = M*patch_properties.stress_ss;
			double eff_oi = M*patch_properties.stress_oi;
			M_Assert(eff_oi>eff_ss, "resolved shear stress on the shearing plane at the onset of instability is lower than at the stable state");
#endif

			if(not continue_shear_test())
			{
				patch_properties.failed = true;

#ifdef DEBUG
				std::cout << "* Warning: probing system failed" << std::endl;
#endif
			}

			data_patch.push_back(patch_properties);

			delete probing_system;

			for(auto &element : patch_elements_copy)
				delete element;

			patch_elements_copy.clear();

			probring_solver.clear();

		} // reference element loop
	} // theta loop
}


/*! Snapshot of patch local properties. */
template<int dim>
class PatchPropertiesSnapshot
{
  public:

	/*! This struct is used for output purposes only. It converts the struct
	 * @ref PatchPropertiesTensorial<dim>, which contains complex objects, into
	 * a plain struct of scalar values. In this way, it can be easily written
	 * into, e.g., hdf5 datasets. */
	struct DataRow
	{
		unsigned int ref_element = 0;
		/*!<  Number of the element used as a reference when creating a patch. It
		 * can be used a patch id to distinguish it from other patches. */

		float theta = 0.;
		/*!< Orientation of the shear test performed on the patch. This angle is
		 * the same used in @ref utils::tensor::schmid. */

		float stress_ss_00 = 0.;
		/*!< Component xx of the stress tensor averaged over the elements composing
		 * the patch. The operation is performed in the ss-state. */

		float stress_ss_11 = 0.;
		/*!< Component yy of the stress tensor averaged over the elements composing
		 * the patch. The operation is performed in the ss-state. */

		float stress_ss_01 = 0.;
		/*!< Component xy of the stress tensor averaged over the elements composing
		 * the patch. The operation is performed in the ss-state. */

		float stress_ee_00 = 0.;
		/*!< Component xx of the stress tensor averaged over the elements composing
		 * the patch. The operation is performed in the ee-state. */

		float stress_ee_11 = 0.;
		/*!< Component yy of the stress tensor averaged over the elements composing
		 * the patch. The operation is performed in the ee-state. */

		float stress_ee_01 = 0.;
		/*!< Component xy of the stress tensor averaged over the elements composing
		 * the patch. The operation is performed in the ee-state. */

		float stress_oi_00 = 0.;
		/*!< Component xx of the stress tensor averaged over the elements composing
		 * the patch. The operation is performed in the oi-state. */

		float stress_oi_11 = 0.;
		/*!< Component yy of the stress tensor averaged over the elements composing
		 * the patch. The operation is performed in the oi-state. */

		float stress_oi_01 = 0.;
		/*!< Component xy of the stress tensor averaged over the elements composing
		 * the patch. The operation is performed in the oi-state. */

		float energy_el_ss = 0.;
		/*!< Elastic energy averaged over the elements composing the patch. The
		 * operation is performed in the ss-state. */

		float energy_el_oi = 0.;
		/*!< Elastic energy averaged over the elements composing the patch. The
		 * operation is performed in the oi-state. */

		float energy_el_ee = 0.;
		/*!< Elastic energy averaged over the elements composing the patch. The
		 * operation is performed in the ee-state. */

		float energy_conf_ss = 0.;
		/*!< Configurational energy averaged over the elements composing the patch. The
		 * operation is performed in the ss-state. */

		float energy_conf_oi = 0.;
		/*!< Configurational energy averaged over the elements composing the patch. The
		 * operation is performed in the oi-state. */

		float energy_conf_ee = 0.;
		/*!< Configurational energy averaged over the elements composing the patch. The
		 * operation is performed in the ee-state. */

		float shear_strain_oi = 0.;
		/*!< Elastic shear strain applied to the patch, with shear orientation
		* @ref theta, to reach the oi-state from the initial ss-state. */

		float x = 0.;
		/*!< x-coordinate of the center of the patch. */

		float y = 0.;
		/*!< y-coordinate of the center of the patch. */

		bool failed = false;
		/*!< Whether the patch has undergone mechanical failure during the shear
		* test. This can be the case if, e.g., strain-softening is present in the
		* system. */

	};

	PatchPropertiesSnapshot(
		const system::System<dim> &system,
		std::string monitor_mag_,
		double desired_target_,
		double recorded_target_,
		unsigned int N_,
		const std::vector<double> &theta_list,
		bool optimized,
		bool do_ee)
		:
		recorded_mag("patches"),
		monitor_name(monitor_mag_),
		desired_target(desired_target_),
		recorded_target(recorded_target_),
		output_index(system.history->index()),
		N(N_)
	{

		/*! Take and store the snapshot from the input system.
		 *
		 * @param system from where the patches are created and analyzed
		 * @param data_patch vector containing the properties measured for each
		 * patch
		 * @param N patch size (see @ref select_patch_elements)
		 * @param theta_list vector with the shear orientations along which shear
		 * tests are to be performed
		 * @param do_ee if false, compute only the ss- and oi-state. If true,
		 * compute also the ee-state
		 */

		std::vector<PatchPropertiesTensorial<dim>> data_;

		if(optimized)
			analyze_patch_ensemble_opt(data_, system, N, theta_list, do_ee);
		else
			analyze_patch_ensemble(data_, system, N, theta_list, do_ee);

		for(auto &d : data_)
		{
			DataRow row;
			row.ref_element = d.ref_element;
			row.theta = d.theta;
			row.stress_ss_00 = d.stress_ss[0][0];
			row.stress_ss_11 = d.stress_ss[1][1];
			row.stress_ss_01 = d.stress_ss[0][1];
			row.stress_oi_00 = d.stress_oi[0][0];
			row.stress_oi_11 = d.stress_oi[1][1];
			row.stress_oi_01 = d.stress_oi[0][1];
			row.shear_strain_oi = d.resolved_elastic_shear_strain_oi;
			row.stress_ee_00 = d.stress_ee[0][0];
			row.stress_ee_11 = d.stress_ee[1][1];
			row.stress_ee_01 = d.stress_ee[0][1];
			row.energy_el_ss = d.energy_el_ss;
			row.energy_el_oi = d.energy_el_oi;
			row.energy_el_ee = d.energy_el_ee;
			row.energy_conf_ss = d.energy_conf_ss;
			row.energy_conf_oi = d.energy_conf_oi;
			row.energy_conf_ee = d.energy_conf_ee;
			row.x = d.coords[0];
			row.y = d.coords[1];
			row.failed = d.failed;

			data.push_back(row);
		}
	};

	std::vector<DataRow> data;
	/*!< Container to store the recorde data. */

	std::string recorded_mag;
	/*!< Name of the field storaged in @ref data. */

	std::string monitor_name;
	/*!< Name of the magnitude used to check whether the snapshot should be taken
	 * or not. */

	double desired_target;
	/*!< Value of the @ref monitor_name at which we desired to take the
	 * snapshot. */

	double recorded_target;
	/*!< Value of the @ref monitor_name at which the snapshot is actually
	 * taken. */

	unsigned int output_index;
	/*!< Global event index from the @ref event::History at which the snapshot
	 * is taken. */

	unsigned int N;
	/*!< Patch size (see @ref select_patch_elements) */

};

} // namespace patches

} // namespace mepls

#endif //__PATCHES_H_
