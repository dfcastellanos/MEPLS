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

#ifndef __SOLVER_H_
#define __SOLVER_H_

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
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/fe/fe_p1nc.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>
#include <string>
#include <map>

#include <mepls/utils.h>
#include <mepls/solver_impl.h>


namespace mepls
{


/*! This namespace contains the elasticity solvers used to compute elastic strain
 * and stress fields induced by the boundary conditions and the plastic (eigen)strain field.
 * The elastic fields are computed by solving the stress equilibrium equation using the Finite
 * Element Method. The solvers are implemented using the <a href="https://www.dealii.org">deal
 * .II</a> library.
 */
namespace elasticity_solver
{

/*! This enum contains values that define the role of the external load. */
enum ControlMode
{
	traction = 0, /*!< The load value controls the applied traction. */
	displacement = 1, /*!< The load value controls the applied displacement. */
	stress = 0, /*!< The load value controls the applied external stress. */
	strain = 1 /*!< The load value controls the applied total strain. */
};

/*!
 * This abstract class provides a common interface for different specialized
 * solver classes.
 *
 * The solver classes compute elastic strain and stress fields induced by
 * the external loading conditions and the plastic strain field (represented by
 * an eigenstrain field). The elastic fields are computed by solving the stress
 * equilibrium equation using the Finite Element Method (FEM) with the
 * <a href="https://www.dealii.org">deal.II</a> library. The solution is
 * computed assuming linear elastic behavior over a regular quadrilateral grid.
 *  Each finite element is associated with a mesoscale element of the system
 *  (see @ref system::System<dim>::elements). The elastic fields are averaged
 *  over each finite element, so a single elastic strain and stress tensor is
 *  associated with each element.
 *
 * @note while dealII denotes each element as a "cell" and represents them by
 * cell objects, the solver classes here described will alternatively refer to
 * them simply as "elements." The elements will be referred to by a 1D index,
 * which is given by the order in which dealII iterates over the existing cells.
 * This index is the same used in the vector of mesoscale elements, @ref
 * system::System<dim>::elements. Moreover, we consider that any property or
 * field defined over the elements is element-wise uniform (as, e.g., elastic
 * stiffness or eigenstrain).
 *
 * The computations are performed using linear shape functions, which increase
 * computational performance at the expense of numerical accuracy in the computation
 * of the elastic fields in comparison with to, e.g., quadratic shape functions.
 * However, the role that elastic fields play in this mesoscale framework allows
 * us to use linear shape functions without a reduction in the quality of the final
 * results that is, in the evolution of plastic activity and microstructural
 * properties of disordered materials. The reason is that the stress
 * fields induced by localized (in space and time) plastic activity is long-ranged
 * (see Eshelby inclusion), and short-range details affecting a small number
 * of elements in the nearest neighborhood of a plastic event are less
 * important than far-field asymptotics. Moreover, the mentioned far-field
 * asymptotics of elastic fields converges to the right solution independently
 * of short-range details. On the other hand, the impact that the loss
 *  of accuracy in the short-range fields might have in the dynamic evolution
 *  of the system is palliated by the fact that we deal with a disordered
 *  material, with a high degree of heterogeneity in its plastic (and possibly
 *  elastic too) properties with stochastic evolution rules. Such disorder has a
 * stronger effect in the noisy plastic activity and response to external loading.
 */
template<int dim>
class Solver
{

	/*!< This struct contains references to the data necessary to assemble FEM
	 * solvers that share the same elastic properties. Its purpose is to easily give
	 * access to existing assembly data to new solvers in order to avoid repeating
	 *  computations. */
	struct DataForAssembling
	{
		DataForAssembling(
			const std::vector<dealii::SymmetricTensor<4, dim>> &C_,
			const std::vector<dealii::FullMatrix<double>> &cell_matrix_assembly_data_,
			const std::vector<std::vector<dealii::Vector<double>>> &unitary_eigenstrains_rhs_)
			:
			C(C_),
			cell_matrix_assembly_data(cell_matrix_assembly_data_),
			unitary_eigenstrains_rhs(unitary_eigenstrains_rhs_)
		{
		};

		const std::vector<dealii::SymmetricTensor<4, dim>> &C;
		const std::vector<dealii::FullMatrix<double>> &cell_matrix_assembly_data;
		const std::vector<std::vector<dealii::Vector<double>>> &unitary_eigenstrains_rhs;
	};


  public:

	Solver(unsigned int Nx, unsigned int Ny)
		:
		Nx(Nx),
		Ny(Ny),
		dof_handler(triangulation),
		fe(dealii::FE_Q<dim>(element_order), dim),
		quadrature_formula(element_order + 1),
		fe_values(fe, quadrature_formula,
				  dealii::update_values | dealii::update_gradients |
				  dealii::update_quadrature_points | dealii::update_JxW_values)
	{
		/*! Constructor.
		  *
		  * @param Nx number of elements in the x-direction
		  * @param Ny number of elements in the y-direction
		  */

		load = 0.;
		already_assembled = false;

		// create a structured quadrilateral mesh
		std::vector<unsigned int> repetitions(2);
		repetitions[0] = Nx;
		repetitions[1] = Ny;
		dealii::GridGenerator::subdivided_hyper_rectangle(triangulation, repetitions,
														  dealii::Point<2>(0.0, 0.0),
														  dealii::Point<2>(Nx * element_length,
																		   Ny * element_length),
														  true);

		dof_handler.distribute_dofs(fe);

		// allocate the containers for the different fields
		solution.reinit(dof_handler.n_dofs());
		local_eigenstrain.resize(triangulation.n_active_cells());
		stress.resize(triangulation.n_active_cells());
		strain.resize(triangulation.n_active_cells());
		elastic_strain.resize(triangulation.n_active_cells());
		deformation_gradient.resize(triangulation.n_active_cells());
		C.resize(triangulation.n_active_cells());

		impl::setup_default_elastic_properties(C);
		impl::average_shape_function_gradients<dim>(*this);
		impl::make_element_maps<dim>(*this);

		const unsigned int dofs_per_cell = fe.dofs_per_cell;
		cell_matrix_assembly_data.resize(triangulation.n_active_cells(),
										 dealii::FullMatrix<double>(dofs_per_cell, dofs_per_cell));
	}

	virtual ~Solver()
	{
		/*! Virtual destructor. */

		dof_handler.clear();
	}

	unsigned int get_n_elements() const
	{
		/*! Get the number of elements in the FEM mesh. */

		return triangulation.n_active_cells();
	}

	void set_elastic_properties(unsigned int n, dealii::SymmetricTensor<4, dim> const &input_C)
	{
		/*! Set the stiffness tensor of the element number `n`. */

		M_Assert(n < get_n_elements(), "");M_Assert(n < C.size(), "");

		set_elastic_properties_impl(n, input_C);
	}


	void set_stress(unsigned int n, const dealii::SymmetricTensor<2, dim> &input_stress)
	{
		/*! Set the stress tensor of the element number `n`.
		 *
		 * @note although stress
		 * tensors are computed by the solver, sometimes the stress tensors are
		 * modified outside the solver, e.g., when a mesoscale element has
		 * a certain pre-stress. With this function, we can inform the solver about
		 * the updated stress tensor. In this way, when the solver compues the
		 * external applied stress, it does so taking into account that field. */

		M_Assert(n < get_n_elements(), "");M_Assert(n < stress.size(), "");

		stress[n] = input_stress;
	}

	const std::vector<dealii::SymmetricTensor<2, dim> > &get_stress()
	{
		/*! Get the stress in all the elements. Each element is characterized
		 * by a single average tensor. */

		return stress;
	}

	const std::vector<dealii::Tensor<2, dim> > &get_deformation_gradient()
	{
		/*! Get the deformation gradient in all the elements. Each element is characterized
		 * by a single average tensor. */

		return deformation_gradient;
	}

	const std::vector<dealii::SymmetricTensor<2, dim> > &get_elastic_strain()
	{
		/*! Get the elastic strain in all the elements. Each element is
		 * characterized by a single average tensor. */

		return elastic_strain;
	}

	const std::vector<dealii::SymmetricTensor<2, dim> > &get_local_eigenstrain()
	{
		/*! Get the local eigenstrain in all the elements. Each element is
		 * characterized by a single tensor, since added eigenstrains are
		 * assumed as uniform over the elements (see @ref add_eigenstrain).
		 *
		 * @note the eigenstrain is generated outside the solver and added to it
		 * using @ref add_eigenstrain, therefore this function only returns
		 * what it has already been given to the solver but does not make
		 * any further computation. */

		return local_eigenstrain;
	}

	DataForAssembling get_data_for_assembling() const
	{
		/*! Get the data used for assembling the FEM solver. This data can
		 * be used to assemble other solvers without recomputing it (see @ref
		 * ShearBoundary<dim>::copy_assembly). */

		return DataForAssembling(C, cell_matrix_assembly_data, unitary_eigenstrains_rhs);
	}

	void add_eigenstrain(unsigned int n, const dealii::SymmetricTensor<2, dim> &eigenstrain)
	{
		/*! Add a local eigenstrain increment that is uniform through the element
		 * number n. The eigenstrain of the rest of the elements is not
		 * changed. */

		M_Assert(already_assembled, "");M_Assert(n < get_n_elements(), "");

		add_eigenstrain_impl(n, eigenstrain);
	}

	void clear()
	{
		/*! Reset the solver to a clean state, where no eigenstrain and no
		 * load is present. The solver remains assembled so it can be used
		 * again. */

		assert(already_assembled);

		load = 0.;
		solution.reinit(dof_handler.n_dofs());
		total_eigenstrain.clear();

		for(unsigned int i = 0; i < triangulation.n_active_cells(); ++i)
		{
			local_eigenstrain[i].clear();
			deformation_gradient[i].clear();
			elastic_strain[i].clear();
			strain[i].clear();
			stress[i].clear();
		}

		clear_impl();
	}

	void solve()
	{
		/*! Solve the elastic equilibrium equation. i.e., compute the elastic
		 * strain and stress fields arising from the external load and the
		 * plastic (eigen)strain field. */

		solve_impl();
	}


	void write_vtu(const std::string &filename)
	{
		/*! Write the displacement, stress and stiffness fields into a VTU
		 * file. */

		std::cout << filename << std::endl;

		impl::Postprocessor<dim> postprocessor(stress, C, cell_to_element);

		dealii::DataOut<dim> data_out;
		data_out.attach_dof_handler(dof_handler);
		data_out.add_data_vector(solution, postprocessor);

		data_out.build_patches();

		std::ofstream output(filename.c_str());
		data_out.write_vtu(output);
	};

	unsigned int get_Nx() const
	{
		/*! Get the number of elements in the x-direction. */

		return Nx;
	};

	unsigned int get_Ny() const
	{
		/*! Get the number of elements in the y-direction. */

		return Ny;
	};

	double get_load() const
	{
		/*! Get the current value of the external load. */

		return load;
	};

	virtual void setup_and_assembly() = 0;
	/*!< Assemble the FEM system. */

	virtual void add_load_increment(double load_increment) = 0;
	/*!< Update the load value by adding the input increment. */

	virtual double get_external_stress() = 0;
	/*!< Get the external stress induced by the load. */

	virtual double get_total_strain() = 0;
	/*!< Get the total strain induced by the load mechanism and the plastic strain field. */

  private:
  	double element_length = 1.;
	/*!< Length of the side of the finite elements. The value is set to 1.0
	 * since it defines the units of length in the framework. */

	unsigned int n_components = 3 * (dim - 1);
	/*!< Number of independent component when writint the symmetric strain and
	 * stress tensors in Voigt's notation. In this ways, the value is 3 if dim=2
	 * and 6 if dim=3. */

	unsigned int element_order = 1;
	/*!< Order of the shape functions. */
	
  protected:

	unsigned int Nx;
	/*!< Number of elements in the x-direction. */

	unsigned int Ny;
	/*!< Number of elements in the y-direction. */

	double load;
	/*!< Value of the external load. The load represents the mechanism externally
	 * driving the system. See @ref mepls::elasticity_solver::ControlMode. */

	bool already_assembled;
	/*!< Bool indicating whether the solver has already been assembled or not. */

	dealii::Vector<double> solution;
	/*!< Displacement field. The values represent the displacement of each degree
	 * of freedom. See <a href="https://www.dealii.org">deal.II</a> documentation
	 * for details. */

	dealii::Triangulation<dim> triangulation;
	/*!< Finite element grid. See <a href="https://www.dealii.org">deal.II</a>
	 * documentation for details. */

	dealii::DoFHandler<dim> dof_handler;
	/*!< Object managing the degrees of freedom. See
	 * <a href="https://www.dealii.org">deal.II</a> documentation for details. */

	dealii::FESystem<dim> fe;
	/*!< Vector-valued finite element object. See
	 * <a href="https://www.dealii.org">deal.II</a> documentation for details. */

	dealii::QGauss<dim> quadrature_formula;
	/*!< The Gauss-Legendre family of quadrature rules for numerical
	 * integration. See <a href="https://www.dealii.org">deal.II</a>
	 * documentation for details. */

	dealii::FEValues<dim> fe_values;
	/*!< Finite element evaluated in the quadrature points of a finite element.
	 * See <a href="https://www.dealii.org">deal.II</a> documentation for
	 * details. */

	std::vector<dealii::FullMatrix<double>> cell_matrix_assembly_data;
	/*!< Vector containing the results of the assembly of each individual element.
	 * See <a href="https://www.dealii.org">deal.II</a> documentation for
	 * details. */

	std::vector<std::vector<double> > av_shape_grads_coeff;
	/*!< Coefficients used to compute the average of the displacement gradient
	 * in an element as a linear combination of the displacements of the degrees
	 * of freedom of the element. See @ref impl::get_av_displacement_gradient for
	 * details . */

	std::vector<typename dealii::DoFHandler<dim>::active_cell_iterator> element_to_cell;
	/*!< Within this solver class, the finite elments are indexed according to an
	 * integer. Each entry of this vector corresponds to that integer, and the
	 * stored cell object corresponds to the dealII's representation of it. */

	std::map<typename dealii::DoFHandler<dim>::active_cell_iterator, unsigned int> cell_to_element;
	/*!< Map between dealII's cells and element number. It is the inverse of
	 * @ref element_to_cell. */

	std::map<unsigned int, std::vector<typename dealii::DoFHandler<dim>::active_cell_iterator> > boundary_map;
	/*!< dealII denotes each domain's boundary by a boundadry id, which is an
	 * integer value. This map relates boundary id's to the cells lying at those
	 * boundaries. */

	std::vector<dealii::SymmetricTensor<4, dim>> C;
	/*!< Vector containing the stiffness tensor of each element. */

	std::vector<dealii::SymmetricTensor<2, dim> > elastic_strain;
	/*!< Vector containing the elastic strain of each element. */

	std::vector<dealii::SymmetricTensor<2, dim> > stress;
	/*!< Vector containing the elastic stress of each element. */

	std::vector<dealii::Tensor<2, dim> > deformation_gradient;
	/*!< Vector containing the deformation gradient of each element. */

	std::vector<dealii::SymmetricTensor<2, dim> > strain;
	/*!< Vector containing the strain of each element. */

	std::vector<dealii::SymmetricTensor<2, dim> > local_eigenstrain;
	/*!< Vector containing the total eigenstrain added to each element. */

	dealii::SymmetricTensor<2, dim> total_eigenstrain;
	/*!< Total eigenstrain added to system. */

	std::vector<std::vector<dealii::Vector<double>>> unitary_eigenstrains_rhs;
	/*!< This structure contains the assembly data necessary to convert any
	 * eigenstrain field into body forces. See
	 * @ref impl::assemble_unitary_eigenstrains and @ref impl::add_eigenstrain.*/

  private:

	virtual void add_eigenstrain_impl(
		unsigned int element, const dealii::SymmetricTensor<dim, 2> &eigenstrain) = 0;
	/*!< This function provides the implementation of
	 * @ref Solver<dim>::add_eigenstrain. */

	virtual void clear_impl() = 0;
	/*!< This function provides the implementation of
	 * @ref Solver<dim>::clear. */

	virtual void solve_impl() = 0;

	/*!< This function provides the implementation of @ref Solver<dim>::solve. */

	virtual void set_elastic_properties_impl(
		unsigned int element, dealii::SymmetricTensor<4, dim> const &CC)
	{
		/* This function provides the implementation of
		 * @ref Solver<dim>::set_elastic_properties. */

		C[element] = CC;
	};

	// declare as friends some functions in the impl namespace, since they are
	// the ones actually implementing the solver and for it they require
	// manipulating its private members
	template<int D>
	friend void impl::calculate_strain(dealii::Vector<double> &solution, Solver<D> &solver);
	template<int D>
	friend void impl::calculate_stress(Solver<D> &solver);
	template<int D>
	friend void impl::assemble_unitary_eigenstrains(Solver<D> &solver);
	template<int D>
	friend void impl::assemble_unitary_eigenstrains_unique(Solver<D> &solver);
	template<int D>
	friend void impl::assemble_eigenstrain_rhs(
		dealii::SymmetricTensor<2, D> &eigenstrain,
		dealii::Vector<double> &cell_rhs,
		typename dealii::DoFHandler<D>::active_cell_iterator &cell,
		Solver<D> &solver);
	template<int D>
	friend void impl::assemble_cell_system_matrix(
		typename dealii::DoFHandler<D>::active_cell_iterator &cell,
		dealii::FEValues<D> &fe_values,
		dealii::FullMatrix<double> &cell_matrix,
		std::vector<dealii::SymmetricTensor<2, D> > &symmetric_gradient,
		dealii::SymmetricTensor<4, D> &C,
		Solver<D> &solver);
	template<int D>
	friend void impl::assemble_system_matrix(
		dealii::SparseMatrix<double> &system_matrix,
		dealii::ConstraintMatrix &constraints,
		dealii::Vector<double> &rhs,
		Solver<D> &solver);
	template<int D>
	friend void impl::assemble_boundary_tractions(
		unsigned int boundary_id,
		dealii::Tensor<1, D> &traction,
		dealii::Vector<double> &rhs,
		Solver<D> &solver);
	template<int D>
	friend void impl::add_eigenstrain(
		const unsigned int &element,
		const dealii::SymmetricTensor<2, D> &eigenstrain,
		dealii::Vector<double> &rhs,
		dealii::ConstraintMatrix &constraints,
		Solver<D> &solver);
	template<int D>
	friend void impl::make_element_maps(Solver<D> &solver);
	template<int D>
	friend void impl::average_shape_function_gradients(Solver<D> &solver);
	template<int D>
	friend void impl::fix_node(
		const std::vector<double> &node_coordinates,
		const unsigned int &dof,
		dealii::ConstraintMatrix &constraints,
		Solver<D> &solver);
};


/*! This solver computes elastic fields with periodic boundary conditions. The
 * external loading conditions are pure shear with principal axes oriented along
 * \f$ \pm \pi/4\f$. The system is driven under strain-controlled conditions.
 *
 *   Since we are dealing with periodic boundary condition, the strain
 *   induced in the system cannot be contolled by prescribing displacement on
 *   the surfaces. As shown by, e.g., @cite MICHEL1999109, the solution can be obtained by
 *   considering that the external load induces an average strain in the system, over which a
 *   periodic fluctuation is added. While the average deformation is obtained
 *    straightforwardly from the loading mode, the fluctuations within the
 *    system induced by the load (as is the case if elastic heterogeneities are
 *    present) are obtained by solving an equivalent problem, namely a
 *    homogeneous eigenstrain field compatible with the imposed average deformation,
 *    and under periodic boundary conditions. Such periodic solution corresponds
 *    to a fluctuation around the average deformation. On the other hand, the
 *    deformation induced by the eigenstrain consequence of plastic deformation
 *    can also be obtained under the constraint of periodic boundary
 *    conditions, and superposed to the deformation induced by the loading. */
template<int dim>
class LeesEdwards: public Solver<dim>
{

	/*! This struct stores the solution computed with an external load value of 1.0
	 * when no plastic (eigen)strain is present. */
	struct SolutionToUnitLoad
	{
		dealii::SymmetricTensor<2, dim> unit_external_strain;
		/*!< Externally applied strain. */

		dealii::Vector<double> displacement;
		/*!< Displacement field, associated to each degree of freedom in the
		 * system. */

		dealii::Vector<double> rhs;
		/*!< Right-hand-side of the linear system of equations created with the
		 * FEM to compute the displacement field. It corresponds to body forces
		 *  associated with each degree of freedom in the system. */

		std::vector<dealii::SymmetricTensor<2, dim> > stress;
		/*!< Vector containing the stress of each element. */

		std::vector<dealii::SymmetricTensor<2, dim> > elastic_strain;
		/*!< Vector containing the elastic strain of each element. */
	};

	/*! This struct contains the information about the internal state of the
	 * solver, so the state of a solver can be modified and then reverted back
	 * to a previous one.
	 *
	 * @note this object stores only those data members that are cleaned when
	 * @ref Solver<dim>::clear() is called. */
	struct State
	{
		double load;
		dealii::Vector<double> solution;
		std::vector<dealii::SymmetricTensor<2, dim> > elastic_strain;
		std::vector<dealii::SymmetricTensor<2, dim> > stress;
		std::vector<dealii::Tensor<2, dim> > deformation_gradient;
		std::vector<dealii::SymmetricTensor<2, dim> > strain;
		std::vector<dealii::SymmetricTensor<2, dim> > local_eigenstrain;
		dealii::SymmetricTensor<2, dim> total_eigenstrain;
		dealii::Vector<double> global_displacement;
		dealii::Vector<double> rhs_load;
		dealii::Vector<double> rhs_eigenstrain;
		dealii::Vector<double> system_rhs;
		dealii::SymmetricTensor<2, dim> external_strain;
		dealii::SymmetricTensor<2, dim> external_stress;
		std::vector<dealii::SymmetricTensor<4, dim>> C;
	};

  public:

	LeesEdwards(unsigned int Nx, unsigned int Ny, ControlMode control_mode = ControlMode::strain);
	/*!< Constructor.
	*
	* @param Nx number of elements in the x-direction
	* @param Ny number of elements in the y-direction
	*/

	State get_state() const;
	/*!< Get the internal state of the solver. This state can be used to revert
	 * its state at a later stage. */

	void set_state(const State &state);
	/*!< Set the internal state of the solver. This state can be used to revert
	 * its state at a later stage. */

	void reassemble_with_new_stiffness(const dealii::SymmetricTensor<4,dim> &input_C);
	/*!< This function assembles the solver using homogeneous elastic properties
	 * given by the input stiffness tensor. Its purpose is the variation of the
	 * global effective elastic properties during a simulatino run. For this
	 * reason the assembly is optimized to leverage the homogeneity of the
	 * properties and does not change the current plastic field and applied load.
	 * The elastic strain and stress are recomputed accordingly using the new elastic
	 * properties.
	 * @warning Use this function only if the elastic properties are homogeneous. */

	void set_control_mode(ControlMode control_mode);
	/*!< Set the control model (stress or strain-controlled). */



	void setup_and_assembly() override;
	void add_load_increment(double load_increment) override;
	double get_external_stress() override;
	double get_total_strain() override;
	using Solver<dim>::add_eigenstrain;
	using Solver<dim>::solve;
	using Solver<dim>::get_elastic_strain;
	using Solver<dim>::get_stress;
	using Solver<dim>::clear;

	dealii::SymmetricTensor<4, dim> get_global_stiffness();
	/*! Compute the effective system-scale stiffness tensor. */

  private:

	std::string active_loading_mode;
	/*!< Name of the loading mode for which the solutions are computed the
	 * solution to three different loading modes to compute the system-scale
	 * stiffness tensor. */

	void solve_for_case(
		const dealii::Vector<double> &rhs, const dealii::Vector<double> &global_displacement);
	/*!< Solve the FEM problem using the external load defined by
	 * @ref active_loading_mode. */

	virtual void clear_impl() override;
	virtual void solve_impl() override;
	virtual void add_eigenstrain_impl(
		unsigned int element, const dealii::SymmetricTensor<2, dim> &eigenstrain) override;

	void compute_global_effective_stiffness();
	/*!< Use the solution to three different loading modes to compute the
	 * system-scale stiffness tensor. */

	void loading_mode_components(unsigned int &i, unsigned int &j) const;
	/*!< Sets the components of the rank-2 tensor associated with the active
	 * loading model. For xy, i=0 and j=1; For xx, i=0 and j=0;
	 * For yy, i=1 and j=1. */

	void assemble_external_strain_load(
		dealii::SymmetricTensor<2, dim> external_strain,
		dealii::Vector<double> &solution_macro_unit_load,
		dealii::Vector<double> &rhs_unit_load);
	/*!< Compute the average displacement solution and body forces induced
	 * by an external unit load. Since we are using linear elasticity, the
	 * solution for any other load value can be obtained by simply multipliying
	 * this unit solution by the desired non-unit load.
	 *
	 * @note this solution corresponds only to the external stress field induced
	 * by the external load. Therefore, it does not include the effects of the
	 * eigenstrain field. Despite the external pure shear conditions, this
	 * solution is in general non-trivial since the system can have elastic
	 * heterogeneities. */

	using Solver<dim>::deformation_gradient;
	using Solver<dim>::strain;
	using Solver<dim>::elastic_strain;
	using Solver<dim>::stress;
	using Solver<dim>::Nx;
	using Solver<dim>::Ny;
	using Solver<dim>::total_eigenstrain;
	using Solver<dim>::local_eigenstrain;
	using Solver<dim>::unitary_eigenstrains_rhs;
	using Solver<dim>::load;
	using Solver<dim>::solution;
	using Solver<dim>::triangulation;
	using Solver<dim>::dof_handler;
	using Solver<dim>::fe;
	using Solver<dim>::quadrature_formula;
	using Solver<dim>::fe_values;
	using Solver<dim>::already_assembled;
	using Solver<dim>::C;
	using Solver<dim>::av_shape_grads_coeff;
	using Solver<dim>::element_to_cell;
	using Solver<dim>::cell_to_element;
	using Solver<dim>::boundary_map;
	using Solver<dim>::cell_matrix_assembly_data;

	dealii::Vector<double> global_displacement;
	/*!< Displacement at each degree of freedom induced by the external loading.
	 * The final, periodic solution corresponds to a fluctuation around that
	 * average. */

	dealii::Vector<double> rhs_load;
	/*!< Force at each degree of freedom induced by the external loading. */

	dealii::Vector<double> rhs_eigenstrain;
	/*!< Force at each degree of freedom induced by the eigenstrain field. */

	dealii::Vector<double> system_rhs;
	/*!< Total force at each degree of freedom induced by the eigenstrain field.
	 * It's the superposition of the forces induced by the eigenstrian field and
	 * the external load. */

	dealii::SparsityPattern sparsity_pattern;
	/*!< Object storing which elements of the @ref system_matrix are nonzero.
	 * See
	 * <a href="https://www.dealii.org">deal.II</a>.*/

	dealii::ConstraintMatrix constraints;
	/*!< Object containning the periodicity constarins and a fixed node
	 * to avoid rigid body motions. See
	 * <a href="https://www.dealii.org">deal.II</a>. */

	dealii::SparseDirectUMFPACK A_direct;
	/*!< FEM direct solver based on UMFPACK. */

	dealii::SparseMatrix<double> system_matrix;
	/*!< FEM stiffness matrix. */

	std::map<std::string, SolutionToUnitLoad> solution_to_loads;
	/*!< Solution for an external unit load under different loading modes (this
	 * includes the @ref global_displacement, @ref rhs_load and the elastic
	 * fields). The available modes are "xy" (pure shear with direction xy);
	 * "xx" (tension with direction x); "yy" (tension with direction y). Since
	 * we are using linear elasticity, the solution for any other load value
	 * can be obtained by simply multipliying this unit solutions by the desired
	 * non-unit load. */

	dealii::SymmetricTensor<4, dim> global_C_eff;
	/*!< System-scale effective stiffness*/

	dealii::SymmetricTensor<2, dim> external_stress;
	/*<! External stress. */

	dealii::SymmetricTensor<2, dim> external_strain;
	/*<! External strain. */

	ControlMode control_mode;
	/*<! Defines the controlled mode (stress or strain). */

};


template<int dim>
LeesEdwards<dim>::LeesEdwards(unsigned int Nx_, unsigned int Ny_, ControlMode control_mode_)
	:
	Solver<dim>(Nx_, Ny_),
	control_mode(control_mode_)
{
	global_displacement.reinit(dof_handler.n_dofs());
	rhs_eigenstrain.reinit(dof_handler.n_dofs());
	rhs_load.reinit(dof_handler.n_dofs());
	system_rhs.reinit(dof_handler.n_dofs());
	active_loading_mode = "xy";
};


template<int dim>
typename LeesEdwards<dim>::State LeesEdwards<dim>::get_state() const
{
	State state;
	state.load = load;
	state.solution = solution;
	state.elastic_strain = elastic_strain;
	state.stress = stress;
	state.deformation_gradient = deformation_gradient;
	state.strain = strain;
	state.local_eigenstrain = local_eigenstrain;
	state.total_eigenstrain = total_eigenstrain;
	state.C = C;
	state.global_displacement = global_displacement;
	state.rhs_load = rhs_load;
	state.rhs_eigenstrain = rhs_eigenstrain;
	state.system_rhs = system_rhs;
	state.external_strain = external_strain;
	state.external_stress = external_stress;

	return state;
}

template<int dim>
void LeesEdwards<dim>::set_state(const LeesEdwards<dim>::State &state)
{
	load = state.load;
	solution = state.solution;
	elastic_strain = state.elastic_strain;
	stress = state.stress;
	deformation_gradient = state.deformation_gradient;
	strain = state.strain;
	local_eigenstrain = state.local_eigenstrain;
	total_eigenstrain = state.total_eigenstrain;
	C = state.C;
	global_displacement = state.global_displacement;
	rhs_load = state.rhs_load;
	rhs_eigenstrain = state.rhs_eigenstrain;
	system_rhs = state.system_rhs;
	external_strain = state.external_strain;
	external_stress = state.external_stress;
}


template<int dim>
void LeesEdwards<dim>::set_control_mode(ControlMode control_mode_)
{
	unsigned int i, j;
	loading_mode_components(i, j);

	if(control_mode == control_mode_)
	{
		return;
	}
	else if(control_mode == ControlMode::strain)
	{
		// displacement is changed by traction
		load = external_stress[i][j];
	}
	else if(control_mode == ControlMode::stress)
	{
		// traction is changed by displacement
		load = external_strain[i][j];
	}
	else
	{
		abort();
	}

	control_mode = control_mode_;
}


template<int dim>
void LeesEdwards<dim>::assemble_external_strain_load(
	dealii::SymmetricTensor<2, dim> external_strain,
	dealii::Vector<double> &solution_macro_unit_load,
	dealii::Vector<double> &rhs_unit_load)
{
	solution_macro_unit_load.reinit(dof_handler.n_dofs());
	rhs_unit_load.reinit(dof_handler.n_dofs());

	// compute the body forces arising in the system due to the application of
	// the external load (event if this load would lead normally to a homogeneous
	// deformation, in this case we consider the general case in which elastic
	// heterogeneoities are present). The external_strain used here corresponds
	// to the deformation of load value 1.0
	typename dealii::DoFHandler<dim>::active_cell_iterator cell;
	unsigned element = 0.;
	for(cell = dof_handler.begin_active(); cell != dof_handler.end(); ++cell, ++element)
		impl::add_eigenstrain(element, -1 * external_strain, rhs_unit_load, constraints, *this);

	// compute the global displacement associated with the average deformation
	// induced by the load
	std::vector<bool> vertex_touched(triangulation.n_vertices(), false);
	for(cell = dof_handler.begin_active(); cell != dof_handler.end(); ++cell)
		for(unsigned int v = 0; v < dealii::GeometryInfo<dim>::vertices_per_cell; ++v)
			if(vertex_touched[cell->vertex_index(v)] == false)
			{
				// iterating over the cells will make us iterate several times each
				// vertex. We check whether we havel already iterated over it to
				// avoid doing it more than once
				vertex_touched[cell->vertex_index(v)] = true;

				// compute the displacement as strain x position
				auto disp = external_strain * cell->vertex(v);
				solution_macro_unit_load[cell->vertex_dof_index(v, 0)] = disp[0];
				solution_macro_unit_load[cell->vertex_dof_index(v, 1)] = disp[1];
			}
}


template<int dim>
void LeesEdwards<dim>::setup_and_assembly()
{
	// create PBC constraints and fix left-bottom node to avoid
	// rigid body motions
	constraints.clear();
	std::vector<double> node = {0., 0.};
	impl::fix_node<dim>(node, 0, constraints, *this);
	impl::fix_node<dim>(node, 1, constraints, *this);
	dealii::DoFTools::make_periodicity_constraints(dof_handler, 0, 1, 0, constraints);
	dealii::DoFTools::make_periodicity_constraints(dof_handler, 2, 3, 1, constraints);
	constraints.close();


	// assemble FEM system
	system_matrix.clear();
	dealii::DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
	dealii::DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, true);
	sparsity_pattern.copy_from(dsp);
	system_matrix.reinit(sparsity_pattern);

	dealii::FEValues<dim> fe_values(fe, quadrature_formula,
									dealii::update_values | dealii::update_gradients |
									dealii::update_quadrature_points | dealii::update_JxW_values);
	const unsigned int dofs_per_cell = fe.dofs_per_cell;
	std::vector<dealii::SymmetricTensor<2, dim> > symmetric_gradient(dofs_per_cell);

	for(unsigned int element = 0; element < triangulation.n_cells(); ++element)
		impl::assemble_cell_system_matrix(element_to_cell[element], fe_values,
										  cell_matrix_assembly_data[element], symmetric_gradient,
										  C[element], *this);

	std::vector<dealii::types::global_dof_index> local_dof_indices(dofs_per_cell);
	typename dealii::DoFHandler<dim>::active_cell_iterator cell;
	unsigned int element = 0;
	for(cell = dof_handler.begin_active(); cell != dof_handler.end(); ++cell, ++element)
	{
		cell->get_dof_indices(local_dof_indices);
		constraints.distribute_local_to_global(cell_matrix_assembly_data[element],
											   local_dof_indices, system_matrix);
	}

	A_direct.initialize(system_matrix);

	impl::assemble_unitary_eigenstrains<dim>(*this);

	// within this call the solution to the unitary load is computed (this is
	// a fundamental step)

	already_assembled = true;

	compute_global_effective_stiffness();

	// after assembling the system, the loading mode is set to pure shear xy
	active_loading_mode = "xy";
}


template<int dim>
void LeesEdwards<dim>::reassemble_with_new_stiffness(const dealii::SymmetricTensor<4,dim> &input_C)
{
	for(unsigned int element = 0; element < triangulation.n_cells(); ++element)
		Solver<dim>::set_elastic_properties(element, input_C);

	/* ------ stiffness matrix assembly ------*/

	if(not already_assembled)
	{
		constraints.clear();
		std::vector<double> node = {0., 0.};
		impl::fix_node<dim>(node, 0, constraints, *this);
		impl::fix_node<dim>(node, 1, constraints, *this);
		dealii::DoFTools::make_periodicity_constraints(dof_handler, 0, 1, 0, constraints);
		dealii::DoFTools::make_periodicity_constraints(dof_handler, 2, 3, 1, constraints);
		constraints.close();

		dealii::DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
		dealii::DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, true);
		sparsity_pattern.copy_from(dsp);
	}

	system_matrix.clear();
	system_matrix.reinit(sparsity_pattern);

	dealii::FEValues<dim> fe_values(fe, quadrature_formula,
									dealii::update_values | dealii::update_gradients |
									dealii::update_quadrature_points | dealii::update_JxW_values);
	const unsigned int dofs_per_cell = fe.dofs_per_cell;
	std::vector<dealii::SymmetricTensor<2, dim> > symmetric_gradient(dofs_per_cell);

	// assembly only one cell of the FEM system, and use it to assemble the rest
	unsigned int element = 0;
	impl::assemble_cell_system_matrix(element_to_cell[element], fe_values,
									  cell_matrix_assembly_data[element], symmetric_gradient,
									  C[element], *this);

	// copy here because it will be used when we copy the assembly into patch solver
	for(unsigned int element = 1; element < triangulation.n_cells(); ++element)
		cell_matrix_assembly_data[element] = cell_matrix_assembly_data[0];

	std::vector<dealii::types::global_dof_index> local_dof_indices(dofs_per_cell);
	typename dealii::DoFHandler<dim>::active_cell_iterator cell;
	for(cell = dof_handler.begin_active(); cell != dof_handler.end(); ++cell, ++element)
	{
		cell->get_dof_indices(local_dof_indices);
		constraints.distribute_local_to_global(cell_matrix_assembly_data[0],
											   local_dof_indices, system_matrix);
	}

	A_direct.initialize(system_matrix);


	/* ------ eigenstrain RHS assembly ------*/

	impl::assemble_unitary_eigenstrains_unique<dim>(*this);


	/* ------ external load RHS assembly ------*/

	dealii::Vector<double> displacement, rhs;
	dealii::SymmetricTensor<2,dim> external_strain__;
	{
		external_strain__[0][1] = 1; // meaning of load value
		assemble_external_strain_load(external_strain__, displacement, rhs);

		SolutionToUnitLoad solution_to_load;
		solution_to_load.displacement = displacement;
		solution_to_load.rhs = rhs;
		solution_to_load.unit_external_strain = external_strain__;

		solve_for_case(rhs, displacement);
		solution_to_load.stress = stress;
		solution_to_load.elastic_strain = elastic_strain;
		solution_to_loads["xy"] = solution_to_load;
	}


	/* ------ set the right solver state ------*/

	// after assembling the system, the loading mode is set to pure shear xy
	active_loading_mode = "xy";

	// we want to change the elastic properties and recompute the elastic fields
	// using the same eigenstrain and external loads that were already present
	auto &rhs_unit_load = solution_to_loads.at(active_loading_mode).rhs;
	auto &global_displacement_unit_load = solution_to_loads.at(active_loading_mode).displacement;
	for(unsigned int i = 0; i < solution.size(); ++i)
	{
		rhs_load[i] = load * rhs_unit_load[i];
		global_displacement[i] = load * global_displacement_unit_load[i];
	}

	rhs_eigenstrain.reinit(rhs_eigenstrain.size());
	for(unsigned int element = 0; element < triangulation.n_cells(); ++element)
		impl::add_eigenstrain<dim>(element, local_eigenstrain[element], rhs_eigenstrain,
		constraints, *this);

	solve();


	already_assembled = true;
}


template<int dim>
void LeesEdwards<dim>::add_eigenstrain_impl(
	unsigned int element, const dealii::SymmetricTensor<2, dim> &eigenstrain)
{
	M_Assert(element < triangulation.n_active_cells(), "");

	impl::add_eigenstrain<dim>(element, eigenstrain, rhs_eigenstrain, constraints, *this);

	local_eigenstrain[element] += eigenstrain;
	total_eigenstrain += eigenstrain;
}


template<int dim>
dealii::SymmetricTensor<4, dim> LeesEdwards<dim>::get_global_stiffness()
{
	return global_C_eff;
};


template<int dim>
void LeesEdwards<dim>::add_load_increment(double load_increment)
{
	M_Assert(already_assembled, "");

	// the solver is implemented using strain-controlled conditions, however after it
	// has been assembled, we can also do stress-controlled since we know the global effective
	// stiffness. To do this, we will compute a strain-equivalent load, which is the one that we
	// will use for computing rhs_load and global_displacement. If we use strain controleld
	// conditions, the strain-equivalent load is simply the load

	load += load_increment;

	// depending on the active loading mode, we must update a certain component
	// of the applied external strain
	unsigned int i, j;
	loading_mode_components(i, j);

	dealii::SymmetricTensor<2,dim> external_strain_incr;
	dealii::SymmetricTensor<2,dim> external_stress_incr;

	if(control_mode == ControlMode::strain)
	{
		external_strain_incr[i][j] = load_increment;
		external_stress_incr = global_C_eff * external_strain_incr;
	}
	else
	{
		external_stress_incr[i][j] = load_increment;
		external_strain_incr = dealii::invert(global_C_eff) * external_stress_incr;
	}

	external_strain += external_strain_incr;
	external_stress += external_stress_incr;
	double eq_strain_load_increment = external_strain_incr[i][j];

	// the rhs and global displacement are updated using the solution to the
	// unit load problem for the active loading mode
	auto &rhs_unit_load = solution_to_loads.at(active_loading_mode).rhs;
	auto &global_displacement_unit_load = solution_to_loads.at(active_loading_mode).displacement;

	for(unsigned int i = 0; i < solution.size(); ++i)
	{
		rhs_load[i] += eq_strain_load_increment * rhs_unit_load[i];
		global_displacement[i] += eq_strain_load_increment * global_displacement_unit_load[i];
	}
}


template<int dim>
void LeesEdwards<dim>::solve_impl()
{
	// the total rhs (body forces) are the superposition of the forces induced
	// by the eigenstrain field (plastic deformation) and the fluctuations
	// induced by the external load (the global average displacement induced by
	// the load will be added directly to the computed displacement field)
	for(unsigned int i = 0; i < system_rhs.size(); ++i)
		system_rhs[i] = rhs_load[i] + rhs_eigenstrain[i];

	solve_for_case(system_rhs, global_displacement);
}


template<int dim>
void LeesEdwards<dim>::solve_for_case(
	const dealii::Vector<double> &rhs, const dealii::Vector<double> &global_displacement)
{
	// compute the PBC solution to the eigenstrain + fluctuation induced by
	// external load
	A_direct.vmult(solution, rhs);
	constraints.distribute(solution);

	// to that solution, add the average displacement induced by the external load
	for(unsigned int i = 0; i < solution.size(); ++i)
		solution[i] += global_displacement[i];

	impl::calculate_strain<dim>(solution, *this);

	if(control_mode == ControlMode::stress)
		for(unsigned int element=0; element < triangulation.n_active_cells(); ++element)
			strain[element] += total_eigenstrain / double(triangulation.n_active_cells());

	impl::calculate_stress<dim>(*this);
}


template<int dim>
void LeesEdwards<dim>::compute_global_effective_stiffness()
{
	// compute the solution to the unitary loads under the 3 different
	// loading modes considered. The solution to the external load is not
	// trivial even if we impose homogeneous loadings since elastic
	// heterogeneities can be present. Use the solutions to those loading
	// modes to compute the global effective stiffness. Those solutions
	// are stored and will be used to solve other problems when calling
	// solve()

	dealii::SymmetricTensor<2,dim> external_strain__;
	dealii::Vector<double> displacement, rhs;
	{
		external_strain__.clear();
		external_strain__[0][0] = 1;
		assemble_external_strain_load(external_strain__, displacement, rhs);

		SolutionToUnitLoad solution_to_load;
		solution_to_load.displacement = displacement;
		solution_to_load.rhs = rhs;
		solution_to_load.unit_external_strain = external_strain__;

		solve_for_case(rhs, displacement);
		solution_to_load.stress = stress;
		solution_to_load.elastic_strain = elastic_strain;
		solution_to_loads["xx"] = solution_to_load;
		clear();
	}
	{
		external_strain__.clear();
		external_strain__[1][1] = 1;
		assemble_external_strain_load(external_strain__, displacement, rhs);

		SolutionToUnitLoad solution_to_load;
		solution_to_load.displacement = displacement;
		solution_to_load.rhs = rhs;
		solution_to_load.unit_external_strain = external_strain__;

		solve_for_case(rhs, displacement);
		solution_to_load.stress = stress;
		solution_to_load.elastic_strain = elastic_strain;
		solution_to_loads["yy"] = solution_to_load;
		clear();
	}
	{
		external_strain__.clear();
		external_strain__[0][1] = 1; // meaning of load value
		assemble_external_strain_load(external_strain__, displacement, rhs);

		SolutionToUnitLoad solution_to_load;
		solution_to_load.displacement = displacement;
		solution_to_load.rhs = rhs;
		solution_to_load.unit_external_strain = external_strain__;

		solve_for_case(rhs, displacement);
		solution_to_load.stress = stress;
		solution_to_load.elastic_strain = elastic_strain;
		solution_to_loads["xy"] = solution_to_load;
		clear();
	}

	auto &stress_a = solution_to_loads["xx"].stress;
	auto &strain_a = solution_to_loads["xx"].elastic_strain;

	auto &stress_b = solution_to_loads["yy"].stress;
	auto &strain_b = solution_to_loads["yy"].elastic_strain;

	auto &stress_c = solution_to_loads["xy"].stress;
	auto &strain_c = solution_to_loads["xy"].elastic_strain;

	// compute averages to obtain global stiffness
	dealii::SymmetricTensor<2, dim> av_stress_a;
	dealii::SymmetricTensor<2, dim> av_stress_b;
	dealii::SymmetricTensor<2, dim> av_stress_c;
	dealii::SymmetricTensor<2, dim> av_strain_a;
	dealii::SymmetricTensor<2, dim> av_strain_b;
	dealii::SymmetricTensor<2, dim> av_strain_c;

	for(unsigned int n = 0; n < triangulation.n_active_cells(); ++n)
	{
		av_stress_a += stress_a[n];
		av_stress_b += stress_b[n];
		av_stress_c += stress_c[n];
		av_strain_a += strain_a[n];
		av_strain_b += strain_b[n];
		av_strain_c += strain_c[n];
	}
	av_stress_a /= double(triangulation.n_active_cells());
	av_stress_b /= double(triangulation.n_active_cells());
	av_stress_c /= double(triangulation.n_active_cells());
	av_strain_a /= double(triangulation.n_active_cells());
	av_strain_b /= double(triangulation.n_active_cells());
	av_strain_c /= double(triangulation.n_active_cells());

	// rhs (stresses)
	dealii::Vector<double> f(9);

	// unkowns (stiffness)
	dealii::Vector<double> u(9);

	// matrix (elastic strains)
	dealii::FullMatrix<double> M(9, 9);

	auto CC = utils::tensor::compute_voigts_stiffness<dim>(av_stress_a, av_stress_b, av_stress_c,
														   av_strain_a, av_strain_b, av_strain_c, f,
														   u, M);

	global_C_eff = utils::tensor::voigt_to_standard_rank4<dim>(CC);
}


template<int dim>
void LeesEdwards<dim>::loading_mode_components(unsigned int &i, unsigned int &j) const
{
	if(active_loading_mode == "xx")
	{
		i = 0;
		j = 0;
	}
	else if(active_loading_mode == "yy")
	{
		i = 1;
		j = 1;
	}
	else if(active_loading_mode == "xy")
	{
		i = 0;
		j = 1;
	}
}


template<int dim>
double LeesEdwards<dim>::get_external_stress()
{
	// since the load imposes homogeneous elastic fields, we can compute the
	// external stress as the average over the system, instead of on the
	// loaded surface. When the system is traction controlled, the external
	//stress is the load, but when it is displacement controlled we must compute it

	unsigned int i, j;
	loading_mode_components(i, j);

	// compute the external stress as the average stress of the relevant
	// stress tensor component (depending on the active loading mode)
	double S = 0.;
	for(unsigned int n = 0; n < stress.size(); ++n)
		S += stress[n][i][j];

	return S / double(stress.size());
};


template<int dim>
double LeesEdwards<dim>::get_total_strain()
{
	// since the load imposes homogeneous elastic fields, we can compute the
	// total strain as the average over the system, instead of on the
	// loaded surface. When the system is displacement controlled, the total
	// strain is the load, but when it is traction controlled we must compute it
	// taking into account the eigenstrain

	unsigned int i, j;
	loading_mode_components(i, j);

	// compute the external stress as the average stress of the relevant
	// stress tensor component (depending on the active loading mode)
	double S = 0.;
	for(unsigned int n = 0; n < strain.size(); ++n)
		S += strain[n][i][j];

	return S / double(strain.size());
};


template<int dim>
void LeesEdwards<dim>::clear_impl()
{
	system_rhs.reinit(dof_handler.n_dofs());
	rhs_eigenstrain.reinit(dof_handler.n_dofs());
	rhs_load.reinit(dof_handler.n_dofs());
	global_displacement.reinit(dof_handler.n_dofs());
	external_strain.clear();
	external_stress.clear();
}


/*! This solver computes elastic fields with non-periodic boundary conditions under
 * an external applied shear strain. The amlitude of the shear strain and its orientation can be
 * controlled by the user
 * @note this solver is designed to operate on patches (see @ref mepls::patches), but can be used
 * on full systems as any other solver, since it has the standard interface defined by @ref
 * mepls::elasticity_solver::Solver<dim>.
  */
template<int dim>
class ShearBoundary: public Solver<dim>
{

	/*! This class provides the values of the displacements that induce a shear
	 * strain along the orientation theta. */
	class ShearBoundaryFunc: public dealii::Function<dim>
	{
	  public:
		ShearBoundaryFunc(double load = 1., double theta = 0.)
			:
			dealii::Function<dim>(2)
		{
			/*! Constructor.
			 *
			 * @param load set the amplitude of the shear strain.
			 * @param theta set the orientation of the shear strain.  */

			set_shear(load, theta);
		}

		void vector_value(const dealii::Point<dim> &p, dealii::Vector<double> &value) const
		{
			/*! Computes the value of the displacement at the spatial point
			 * p, given the applied strain. */

			value[0] = strain[0][0] * p[0] + strain[0][1] * p[1];
			value[1] = strain[1][0] * p[0] + strain[1][1] * p[1];
		}

		void set_shear(double load_, double theta_)
		{
			/*! Create the strain tensor from its amplitude @ref load and
			 * shear orientation @ref theta. */

			load = load_;
			theta = theta_;
			strain = load * 2. * utils::tensor::make_schmid<dim>(theta);
			// e.g., load*2*M(theta=0) corresponds to strain = [[0,load],[load,0]],
			// thus load=eps_xy
		}

	  private:
		dealii::Tensor<2, dim> strain;
		/*!< Applied shear strain. */

		double theta;
		/*!< Orientation of applied the shear strain. */

		double load;
		/*!< Amplitude of the shear strain. It's not \f$ \gamma \f$ but
		 * \f$ \epsilon_{\rm xy} \f$. */
	};


  public:

	ShearBoundary(unsigned int Nx, unsigned int Ny);
	/*!< Constructor.
	*
	* @param Nx number of elements in the x-direction
	* @param Ny number of elements in the y-direction
	*/

	void setup_and_assembly();

	void copy_assembly(
		const std::vector<unsigned int> &elements_to_assemble, const Solver<dim> &input_solver);
	/*!< Replace the result of the assembly of this solver with the result of
	 * the assemble of the input solver in this way, we can assemble a bigger
	 * system with heterogeneous properties, and then "extract them" and use
	 * them in this solver.
	 *
	 * @note only the system matrix must be completely re-assembled from the
	 * copied data since cells that were not at boundaries in the input solver
	 * now can be at the boundaries (in other words, the different constraints
	 * lead to stiffness matrices with different forms). However, the
	 * eigenstrain assembly remains the same (since the elastic properties are
	 * copied to this solver, and the different possible constraints are
	 * applied to the rhs not during assembly but when adding the rhs
	 * contributions, later that when the assembly takes place). */

	void add_load_increment(double load_increment);
	double get_external_stress();
	double get_total_strain();

	double theta;
	/*!< Orientation of applied the shear strain. */

	dealii::Vector<double> system_rhs;
	/*!< Total force at each degree of freedom induced by the eigenstrain field.
	 * It's the superposition of the forces induced by the eigenstrian field and
	 * the external load. */

	dealii::Vector<double> rhs_unit_load;
	/*!< Force at each degree of freedom induced by the external loading. */

	dealii::SparsityPattern sparsity_pattern;
	/*!< Object storing which elements of the @ref system_matrix are nonzero.
	 * See <a href="https://www.dealii.org">deal.II</a>.*/

	dealii::SparseDirectUMFPACK A_direct;
	/*!< FEM direct solver based on UMFPACK. */

	dealii::SparseMatrix<double> system_matrix;
	/*!< FEM stiffness matrix. */

	dealii::ConstraintMatrix constraints;
	/*!< Object containning the periodicity constarins and a fixed node
	 * to avoid rigid body motions. See
	 * <a href="https://www.dealii.org">deal.II</a>. */

	using Solver<dim>::add_eigenstrain;
	using Solver<dim>::solve;
	using Solver<dim>::clear;
	using Solver<dim>::deformation_gradient;
	using Solver<dim>::strain;
	using Solver<dim>::elastic_strain;
	using Solver<dim>::stress;
	using Solver<dim>::Nx;
	using Solver<dim>::Ny;
	using Solver<dim>::total_eigenstrain;
	using Solver<dim>::local_eigenstrain;
	using Solver<dim>::unitary_eigenstrains_rhs;
	using Solver<dim>::load;
	using Solver<dim>::solution;
	using Solver<dim>::triangulation;
	using Solver<dim>::dof_handler;
	using Solver<dim>::fe;
	using Solver<dim>::quadrature_formula;
	using Solver<dim>::fe_values;
	using Solver<dim>::already_assembled;
	using Solver<dim>::C;
	using Solver<dim>::av_shape_grads_coeff;
	using Solver<dim>::element_to_cell;
	using Solver<dim>::cell_to_element;
	using Solver<dim>::boundary_map;
	using Solver<dim>::cell_matrix_assembly_data;

  private:
	ShearBoundaryFunc shear_boundary;
	/*!< Function used for creating the shear strain boundary conditions. */

	void set_value_shear_boundary_constraint(double load, double theta);
	/*!< Establish the amplitude (load) and orientation of the shear boundary
	 * conditions. */

	void clear_impl();
	void solve_impl();
	void add_eigenstrain_impl(
		unsigned int element, const dealii::SymmetricTensor<2, dim> &eigenstrain);

	void perform_matrix_assembly_with_existing_data(std::vector<dealii::FullMatrix<double>> &cell_matrix_assembly_data);
	/*!< Assemble the FEM stiffness matix @ref system_matrix using the input
	 * assembly data. */
};


template<int dim>
ShearBoundary<dim>::ShearBoundary(unsigned int Nx_, unsigned int Ny_)
	:
	Solver<dim>(Nx_, Ny_)
{
	theta = 0.;

	system_rhs.reinit(dof_handler.n_dofs());
	rhs_unit_load.reinit(dof_handler.n_dofs());
}


template<int dim>
void ShearBoundary<dim>::set_value_shear_boundary_constraint(double load, double theta)
{
	shear_boundary.set_shear(load, theta);

	constraints.clear();
	dealii::VectorTools::interpolate_boundary_values(dof_handler, 0, shear_boundary, constraints);
	dealii::VectorTools::interpolate_boundary_values(dof_handler, 1, shear_boundary, constraints);
	dealii::VectorTools::interpolate_boundary_values(dof_handler, 2, shear_boundary, constraints);
	dealii::VectorTools::interpolate_boundary_values(dof_handler, 3, shear_boundary, constraints);
	constraints.close();
}


template<int dim>
void ShearBoundary<dim>::setup_and_assembly()
{
	// set_value_shear_boundary_constraint must be call before creating the
	// sparsity pattern
	set_value_shear_boundary_constraint(1, theta);

	sparsity_pattern.reinit(dof_handler.n_dofs(), dof_handler.n_dofs(),
							dof_handler.max_couplings_between_dofs());
	dealii::DoFTools::make_sparsity_pattern(dof_handler, sparsity_pattern);
	constraints.condense(sparsity_pattern);
	sparsity_pattern.compress();


	impl::assemble_unitary_eigenstrains<dim>(*this);

	const unsigned int dofs_per_cell = fe.dofs_per_cell;
	std::vector<dealii::SymmetricTensor<2, dim> > symmetric_gradient(dofs_per_cell);
	dealii::FEValues<dim> fe_values(fe, quadrature_formula,
									dealii::update_values | dealii::update_gradients |
									dealii::update_quadrature_points | dealii::update_JxW_values);

	for(unsigned int element = 0; element < triangulation.n_cells(); ++element)
		impl::assemble_cell_system_matrix(element_to_cell[element], fe_values,
										  cell_matrix_assembly_data[element], symmetric_gradient,
										  C[element], *this);


	perform_matrix_assembly_with_existing_data(cell_matrix_assembly_data);

	already_assembled = true;
}


template<int dim>
void ShearBoundary<dim>::perform_matrix_assembly_with_existing_data(
	std::vector<dealii::FullMatrix<double>> &cell_matrix_assembly_data)
{
	M_Assert(cell_matrix_assembly_data.size() == triangulation.n_active_cells(), "");

	rhs_unit_load.reinit(dof_handler.n_dofs());
	system_matrix.reinit(sparsity_pattern);

	const unsigned int dofs_per_cell = fe.dofs_per_cell;
	dealii::Vector<double> fake_cell_rhs(dofs_per_cell);
	fake_cell_rhs = 0;
	dealii::Vector<double> fake_rhs;
	fake_rhs.reinit(dof_handler.n_dofs());
	std::vector<dealii::types::global_dof_index> local_dof_indices(dofs_per_cell);
	typename dealii::DoFHandler<dim>::active_cell_iterator cell;

	// unitary load during assembly, to get the rhs_unit_load vector
	set_value_shear_boundary_constraint(1, theta);

	for(unsigned int n = 0; n < triangulation.n_active_cells(); ++n)
	{
		cell = element_to_cell[n];
		cell->get_dof_indices(local_dof_indices);
		constraints.distribute_local_to_global(cell_matrix_assembly_data[n], fake_cell_rhs,
											   local_dof_indices, system_matrix, rhs_unit_load);
	}

	A_direct.initialize(system_matrix);

	set_value_shear_boundary_constraint(0., theta); // applied load is set to 0
}


template<int dim>
void ShearBoundary<dim>::copy_assembly(
	const std::vector<unsigned int> &input_elements_to_copy_from, const Solver<dim> &input_solver)
{
	M_Assert(already_assembled, "Normal assembly must be performed at least once");

	auto data_for_asesmbly = input_solver.get_data_for_assembling();
	auto &input_cell_matrix_assembly = data_for_asesmbly.cell_matrix_assembly_data;
	auto &input_unitary_eigenstrains_rhs = data_for_asesmbly.unitary_eigenstrains_rhs;
	auto &input_C = data_for_asesmbly.C;

	for(unsigned int n = 0; n < input_elements_to_copy_from.size(); ++n)
	{
		unsigned int e_from = input_elements_to_copy_from[n];
		unsigned int e_to = n;

		cell_matrix_assembly_data[e_to] = input_cell_matrix_assembly[e_from];
		unitary_eigenstrains_rhs[e_to] = input_unitary_eigenstrains_rhs[e_from];
		C[e_to] = input_C[e_from];
	}

	perform_matrix_assembly_with_existing_data(cell_matrix_assembly_data);
}


template<int dim>
void ShearBoundary<dim>::add_eigenstrain_impl(
	unsigned int element, const dealii::SymmetricTensor<2, dim> &eigenstrain)
{
	M_Assert(element < triangulation.n_active_cells(), "");

	impl::add_eigenstrain<dim>(element, eigenstrain, system_rhs, constraints, *this);

	local_eigenstrain[element] += eigenstrain;
	total_eigenstrain += eigenstrain;
}


template<int dim>
void ShearBoundary<dim>::add_load_increment(double load_increment)
{
	M_Assert(already_assembled, "");

	load += load_increment;

	// update the body forces for the new load
	for(unsigned int i = 0; i < system_rhs.size(); ++i)
		system_rhs[i] += load_increment * rhs_unit_load[i];

	// updted the displacements induced by the new load
	set_value_shear_boundary_constraint(load, theta);
}

template<int dim>
double ShearBoundary<dim>::get_external_stress()
{
	// the external stress is here the shear stress on the shear plane,
	// defined by theta

	dealii::SymmetricTensor<2, dim> ext_stress;
	for(auto &tensor : stress)
		ext_stress += tensor;
	ext_stress /= double(stress.size());

	double resolved_shear_stress = ext_stress * utils::tensor::make_schmid<dim>(theta);

	return resolved_shear_stress;
}


template<int dim>
double ShearBoundary<dim>::get_total_strain()
{
	return load;
}


template<int dim>
void ShearBoundary<dim>::solve_impl()
{
	M_Assert(already_assembled, "");

	A_direct.vmult(solution, system_rhs);

	constraints.distribute(solution);

	impl::calculate_strain<dim>(solution, *this);

	impl::calculate_stress<dim>(*this);
}


template<int dim>
void ShearBoundary<dim>::clear_impl()
{
	// reset the displacements to 0
	set_value_shear_boundary_constraint(0., theta);

	// reset the body forces 0
	system_rhs.reinit(dof_handler.n_dofs());
}


} // namespace elasticity_solver
} // namespace mepls

#endif //__SOLVER_H_
