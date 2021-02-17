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

#ifndef __SolverImpl__
#define __SolverImpl__

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/affine_constraints.h>
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


namespace mepls
{
namespace elasticity_solver
{

template<int dim>
class Solver;

enum ControlMode
{
	traction = 0, displacement = 1
};


namespace impl
{

template<int dim>
class Postprocessor: public dealii::DataPostprocessor<dim>
{
  public:

	Postprocessor(
		std::vector<dealii::SymmetricTensor<2, dim> > &cell_stress_,
		std::vector<dealii::SymmetricTensor<4, dim>> &C_,
		std::map<typename dealii::DoFHandler<dim>::active_cell_iterator, unsigned int> &cell_to_element_)
		:
		cell_stress(cell_stress_),
		cell_to_element(cell_to_element_),
		C(C_)
	{
	}

	virtual void evaluate_vector_field(
		const dealii::DataPostprocessorInputs::Vector<dim> &inputs,
		std::vector<dealii::Vector<double> > &computed_quantities) const
	{
		typename dealii::DoFHandler<dim>::active_cell_iterator cell = inputs.template get_cell<dealii::DoFHandler<dim>>();
		unsigned int element = cell_to_element[cell];
		const unsigned int n_quadrature_points = inputs.solution_values.size();

		std::vector<double> G;
		for(auto &c : C)
			G.push_back(c[0][1][0][1]);

		for(unsigned int q = 0; q < n_quadrature_points; ++q)
		{
			for(unsigned int d = 0; d < dim; ++d)
				computed_quantities[q](d) = inputs.solution_values[q](d);

			//				strain_tensor[0][0] = inputs.solution_gradients[q][0][0];
			//				strain_tensor[1][1] = inputs.solution_gradients[q][1][1];
			//				strain_tensor[0][1] = 0.5*(inputs.solution_gradients[q][0][1]+inputs.solution_gradients[q][1][0]);
			//				el_strain_tensor = strain_tensor-local_eigenstrain[element];
			//				stress_tensor = C*el_strain_tensor;

			computed_quantities[q](dim) = cell_stress[element][0][0];
			computed_quantities[q](dim + 1) = cell_stress[element][1][1];
			computed_quantities[q](dim + 2) = cell_stress[element][0][1];
			computed_quantities[q](dim + 3) = G[element];
		}
	}


	virtual std::vector<std::string> get_names() const
	{
		std::vector<std::string> solution_names(dim, "displacement");
		solution_names.push_back("stress_00");
		solution_names.push_back("stress_11");
		solution_names.push_back("stress_01");
		solution_names.push_back("G");

		return solution_names;
	}

	virtual std::vector<dealii::DataComponentInterpretation::DataComponentInterpretation> get_data_component_interpretation() const
	{
		std::vector<dealii::DataComponentInterpretation::DataComponentInterpretation> interpretation(
			dim, dealii::DataComponentInterpretation::component_is_part_of_vector);
		interpretation.push_back(dealii::DataComponentInterpretation::component_is_scalar);
		interpretation.push_back(dealii::DataComponentInterpretation::component_is_scalar);
		interpretation.push_back(dealii::DataComponentInterpretation::component_is_scalar);
		interpretation.push_back(dealii::DataComponentInterpretation::component_is_scalar);
		return interpretation;
	}

	virtual dealii::UpdateFlags get_needed_update_flags() const
	{
		return dealii::update_values | dealii::update_gradients | dealii::update_quadrature_points;
	}

	std::vector<dealii::SymmetricTensor<2, dim> > &cell_stress;
	std::map<typename dealii::DoFHandler<dim>::active_cell_iterator, unsigned int> &cell_to_element;
	std::vector<dealii::SymmetricTensor<4, dim>> &C;
};


template<int dim>
void fix_node(
	const std::vector<double> &node_coordinates,
	const unsigned int &dof,
	dealii::AffineConstraints<double> &constraints,
	Solver<dim> &solver)
{
	dealii::Point<dim> node;

	for(unsigned int i = 0; i < dim; ++i)
		node[i] = node_coordinates[i];

	typename dealii::DoFHandler<dim>::active_cell_iterator cell;

	for(cell = solver.dof_handler.begin_active(); cell != solver.dof_handler.end(); ++cell)
		for(unsigned int v = 0; v < dealii::GeometryInfo<dim>::vertices_per_cell; ++v)
			if(std::sqrt((cell->vertex(v) - node).norm_square()) < 1e-8)
				constraints.add_line(cell->vertex_dof_index(v, dof));
}


template<int dim>
void assemble_boundary_tractions(
	unsigned int boundary_id,
	dealii::Tensor<1, dim> &traction,
	dealii::Vector<double> &rhs,
	Solver<dim> &solver)
{
	const unsigned int dofs_per_cell = solver.fe.dofs_per_cell;
	std::vector<dealii::types::global_dof_index> local_dof_indices(dofs_per_cell);

	const dealii::FEValuesExtractors::Vector displacements(0);
	dealii::Vector<double> cell_rhs(dofs_per_cell);

	// For a given boundary id, assemble the rhs traction contribution of the face of the reference cell that is at that boundary.
	// It is necessary that the boundary id's are the same as the numbering of the faces of the reference cell.
	unsigned int face_id = boundary_id;
	const unsigned int n_faces = dealii::GeometryInfo<dim>::faces_per_cell;
	const dealii::UpdateFlags uf_face(
		dealii::update_values | dealii::update_normal_vectors | dealii::update_JxW_values);
	const dealii::QGauss<dim - 1> qf_face(solver.element_order + 1);
	dealii::FEFaceValues<dim> fe_values_f(solver.fe, qf_face, uf_face);
	const unsigned int n_q_points_f(qf_face.size());

	typename dealii::DoFHandler<dim>::active_cell_iterator cell = solver.dof_handler.begin_active(), endc = solver.dof_handler.end();

	for(; cell != endc; ++cell)
		if(cell->face(face_id)->at_boundary())
		{
			fe_values_f.reinit(cell, face_id);
			cell_rhs.reinit(dofs_per_cell);

			for(unsigned int f_q_point = 0; f_q_point < n_q_points_f; ++f_q_point)
			{
				//	    const dealii::Tensor<1,dim> traction = load * fe_values_f.normal_vector(f_q_point); // for non-square meshes we can get automatically the normal vectors to the surface
				for(unsigned int i = 0; i < dofs_per_cell; ++i)
					cell_rhs(i) += (fe_values_f[displacements].value(i,
																	 f_q_point) * traction) * fe_values_f.JxW(
						f_q_point);
			}

			cell->get_dof_indices(local_dof_indices);
			for(unsigned int i = 0; i < dofs_per_cell; ++i)
				rhs(local_dof_indices[i]) += cell_rhs[i];
		}
}


template<int dim>
void assemble_eigenstrain_rhs(
	dealii::SymmetricTensor<2, dim> &eigenstrain,
	dealii::Vector<double> &cell_rhs,
	typename dealii::DoFHandler<dim>::active_cell_iterator &cell,
	Solver<dim> &solver)
{
	const unsigned int dofs_per_cell = solver.fe.dofs_per_cell;
	const unsigned int n_q_points = solver.quadrature_formula.size();
	solver.fe_values.reinit(cell);
	cell_rhs.reinit(dofs_per_cell);
	const dealii::FEValuesExtractors::Vector displacements(0);
	dealii::SymmetricTensor<2, dim> symmetric_gradient;
	unsigned int element = solver.cell_to_element.at(cell);

	for(unsigned int q = 0; q < n_q_points; ++q)
		for(unsigned int dof = 0; dof < dofs_per_cell; ++dof)
		{
			symmetric_gradient = solver.fe_values[displacements].symmetric_gradient(dof, q);
			cell_rhs(
				dof) += symmetric_gradient * solver.C[element] * eigenstrain * solver.fe_values.JxW(
				q); // eigenstress is constant in each element, then does not depend on the quadrature points
		}
}


template<int dim>
void assemble_unitary_eigenstrains(Solver<dim> &solver)
{
	/* Assemble the "canonical basis" of eigenstrains.
	 * In linear elasticity, any other eigenstrain RHS
	 * contribution will be a linear combinations of these.
	 */

	dealii::Vector<double> cell_rhs;
	dealii::SymmetricTensor<2, dim> eigenstrain;
	std::vector<dealii::Vector<double>> eigenstrain_rhs_at_element(solver.n_components);
	typename dealii::DoFHandler<dim>::active_cell_iterator cell = solver.dof_handler.begin_active(), endc = solver.dof_handler.end();

	unsigned int element = 0;
	for(; cell != endc; ++cell, ++element)
	{
		eigenstrain[0][0] = 1.;
		eigenstrain[1][1] = 0.;
		eigenstrain[0][1] = 0.;
		assemble_eigenstrain_rhs(eigenstrain, cell_rhs, cell, solver);
		eigenstrain_rhs_at_element[0] = cell_rhs;

		eigenstrain[0][0] = 0.;
		eigenstrain[1][1] = 1.;
		eigenstrain[0][1] = 0.;
		assemble_eigenstrain_rhs(eigenstrain, cell_rhs, cell, solver);
		eigenstrain_rhs_at_element[1] = cell_rhs;

		eigenstrain[0][0] = 0.;
		eigenstrain[1][1] = 0.;
		eigenstrain[0][1] = 1.;
		assemble_eigenstrain_rhs(eigenstrain, cell_rhs, cell, solver);
		eigenstrain_rhs_at_element[2] = cell_rhs;

		solver.unitary_eigenstrains_rhs.push_back(eigenstrain_rhs_at_element);
	}
}


template<int dim>
void add_eigenstrain(
	const unsigned int &element,
	const dealii::SymmetricTensor<2, dim> &eigenstrain,
	dealii::Vector<double> &rhs,
	dealii::AffineConstraints<double> &constraints,
	Solver<dim> &solver)
{
	const unsigned int dofs_per_cell = solver.fe.dofs_per_cell;
	std::vector<dealii::types::global_dof_index> local_dof_indices(dofs_per_cell);

	typename dealii::DoFHandler<dim>::active_cell_iterator cell = solver.element_to_cell[element];
	cell->get_dof_indices(local_dof_indices);

	dealii::Vector<double> cell_rhs;
	cell_rhs.reinit(dofs_per_cell);

	for(unsigned int i = 0; i < dofs_per_cell; ++i)
		cell_rhs(i) += solver.unitary_eigenstrains_rhs[element][0](
			i) * eigenstrain[0][0] + solver.unitary_eigenstrains_rhs[element][1](
			i) * eigenstrain[1][1] + solver.unitary_eigenstrains_rhs[element][2](
			i) * eigenstrain[0][1];

	constraints.distribute_local_to_global(cell_rhs, local_dof_indices, rhs);
}


template<int dim>
void make_element_maps(Solver<dim> &solver)
{
	unsigned int n_active_cells = solver.triangulation.n_active_cells();
	typename dealii::DoFHandler<dim>::active_cell_iterator cell = solver.dof_handler.begin_active();
	unsigned int faces_per_cell = dealii::GeometryInfo<dim>::faces_per_cell;

	solver.element_to_cell.resize(n_active_cells);

	for(unsigned int i = 0; i < n_active_cells; ++i)
	{
		solver.element_to_cell[i] = cell;
		solver.cell_to_element[cell] = i;

		for(unsigned int id = 0; id < faces_per_cell; ++id)
			if(cell->face(id)->at_boundary())
				solver.boundary_map[id].push_back(cell);

		++cell;
	}
}


template<int dim>
void average_shape_function_gradients(Solver<dim> &solver)
{
	const unsigned int vertices_per_cell = dealii::GeometryInfo<dim>::vertices_per_cell;
	const unsigned int dofs_per_vertex = 2;

	solver.av_shape_grads_coeff.resize(dofs_per_vertex);
	for(unsigned int i = 0; i < dofs_per_vertex; ++i)
		solver.av_shape_grads_coeff[i].resize(vertices_per_cell, 0.0);

	double Lx = solver.dof_handler.begin_active()->extent_in_direction(0);
	double Ly = solver.dof_handler.begin_active()->extent_in_direction(1);

	solver.av_shape_grads_coeff[0][0] = -0.5 / Lx;
	solver.av_shape_grads_coeff[0][1] = 0.5 / Lx;
	solver.av_shape_grads_coeff[0][2] = -0.5 / Lx;
	solver.av_shape_grads_coeff[0][3] = 0.5 / Lx;
	solver.av_shape_grads_coeff[1][0] = -0.5 / Ly;
	solver.av_shape_grads_coeff[1][1] = -0.5 / Ly;
	solver.av_shape_grads_coeff[1][2] = 0.5 / Ly;
	solver.av_shape_grads_coeff[1][3] = 0.5 / Ly;
}


template<int dim>
void assemble_cell_system_matrix(
	typename dealii::DoFHandler<dim>::active_cell_iterator &cell,
	dealii::FEValues<dim> &fe_values,
	dealii::FullMatrix<double> &cell_matrix,
	std::vector<dealii::SymmetricTensor<2, dim> > &symmetric_gradient,
	dealii::SymmetricTensor<4, dim> &C,
	Solver<dim> &solver)
{
	const unsigned int n_q_points = solver.quadrature_formula.size();
	const unsigned int dofs_per_cell = solver.fe.dofs_per_cell;
	const dealii::FEValuesExtractors::Vector displacements(0);

	cell_matrix = 0;
	fe_values.reinit(cell);

	for(unsigned int q_point = 0; q_point < n_q_points; ++q_point)
	{
		for(unsigned int k = 0; k < dofs_per_cell; ++k)
			symmetric_gradient[k] = fe_values[displacements].symmetric_gradient(k, q_point);

		for(unsigned int i = 0; i < dofs_per_cell; ++i)
			for(unsigned int j = 0; j < dofs_per_cell; ++j)
				cell_matrix(i,
							j) += (symmetric_gradient[i] * C * symmetric_gradient[j]) * fe_values.JxW(
					q_point);
	}
}


template<int dim>
void assemble_system_matrix(
	dealii::SparseMatrix<double> &system_matrix,
	dealii::AffineConstraints<double> &constraints,
	dealii::Vector<double> &rhs,
	Solver<dim> &solver)
{
	unsigned int n_active_cells = solver.triangulation.n_active_cells();
	const unsigned int dofs_per_cell = solver.fe.dofs_per_cell;

	dealii::Vector<double> fake_cell_rhs(dofs_per_cell);
	fake_cell_rhs = 0;

	dealii::Vector<double> fake_rhs;
	fake_rhs.reinit(solver.dof_handler.n_dofs());

	std::vector<dealii::SymmetricTensor<2, dim> > symmetric_gradient(dofs_per_cell);
	dealii::FEValues<dim> fe_values(solver.fe, solver.quadrature_formula,
									dealii::update_values | dealii::update_gradients | dealii::update_quadrature_points | dealii::update_JxW_values);
	std::vector<dealii::types::global_dof_index> local_dof_indices(dofs_per_cell);
	dealii::FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
	typename dealii::DoFHandler<dim>::active_cell_iterator cell;

	for(unsigned int element = 0; element < n_active_cells; ++element)
	{
		cell = solver.element_to_cell[element];

		assemble_cell_system_matrix(cell, fe_values, cell_matrix, symmetric_gradient,
									solver.C[element], solver);

		cell->get_dof_indices(local_dof_indices);

		for(unsigned int i = 0; i < dofs_per_cell; ++i)
			for(unsigned int j = 0; j < dofs_per_cell; ++j)
				system_matrix.add(local_dof_indices[i], local_dof_indices[j], cell_matrix(i, j));

		constraints.condense(system_matrix, rhs);
	}
}


template<int dim>
void get_av_displacement_gradient(
	dealii::Tensor<2, dim> &J,
	std::vector<std::vector<double> > &av_shape_grads_coeff,
	dealii::Vector<double> &solution,
	typename dealii::DoFHandler<dim>::active_cell_iterator &cell)
{
	const unsigned int vertices_per_cell = dealii::GeometryInfo<dim>::vertices_per_cell;
	J.clear();
	for(unsigned int k = 0; k < vertices_per_cell; ++k)
	{
		J[0][0] += av_shape_grads_coeff[0][k] * solution(cell->vertex_dof_index(k, 0));
		J[1][1] += av_shape_grads_coeff[1][k] * solution(cell->vertex_dof_index(k, 1));
		J[0][1] += av_shape_grads_coeff[1][k] * solution(cell->vertex_dof_index(k, 0));
		J[1][0] += av_shape_grads_coeff[0][k] * solution(cell->vertex_dof_index(k, 1));
	}
}


template<int dim>
void get_av_displacement_gradient_dealii_builtin(
	dealii::Tensor<2, dim> &J,
	Solver<dim> &solver,
	typename dealii::DoFHandler<dim>::active_cell_iterator &cell)
{
	// get grandients from dealII (useful if we have higher order shape functions)


	// these objects are created for each cell. They could be reused to increase the performance
	const double n_q_points = solver.quadrature_formula.size();
	std::vector<dealii::Tensor<2, dim>> dJ(n_q_points);
	const dealii::FEValuesExtractors::Vector displacements(0);

	solver.fe_values.reinit(cell);
	solver.fe_values[displacements].get_function_gradients(solver.solution, dJ);
	J.clear();
	for(unsigned int q = 0; q < n_q_points; ++q)
		J += dJ[q];
	J /= double(n_q_points);
}


template<int dim>
void calculate_strain(dealii::Vector<double> &solution, Solver<dim> &solver)
{
	typename dealii::DoFHandler<dim>::active_cell_iterator cell = solver.dof_handler.begin_active(), endc = solver.dof_handler.end();
	dealii::Tensor<2, dim> J;
	dealii::Tensor<2, dim> J_sym;
	dealii::Tensor<2, dim> identity;
	identity[0][0] = 1;
	identity[1][1] = 1;

	unsigned int element = 0;

	for(; cell != endc; ++cell, ++element)
	{
		get_av_displacement_gradient(J, solver.av_shape_grads_coeff, solution, cell);
		//		get_av_displacement_gradient_dealii_builtin(J, solver, cell);

		J_sym = J + transpose(J);

		solver.strain[element][0][0] = 0.5 * J_sym[0][0];
		solver.strain[element][1][1] = 0.5 * J_sym[1][1];
		solver.strain[element][0][1] = 0.5 * J_sym[0][1];

		solver.deformation_gradient[element] = identity + J;
	}
}


template<int dim>
void calculate_stress(Solver<dim> &solver)
{
	for(unsigned int element = 0; element < solver.triangulation.n_active_cells(); ++element)
	{
		solver.elastic_strain[element] = solver.strain[element] - solver.local_eigenstrain[element];
		solver.stress[element] = solver.C[element] * solver.elastic_strain[element];
	}
}


template<int dim>
void setup_default_elastic_properties(std::vector<dealii::SymmetricTensor<4, dim>> &C)
{
	dealii::SymmetricTensor<2, 3> CC;
	double nu = 0.3;
	double G = 1.;
	double E = 2 * G * (1 + nu);
	double c = E / (1. - nu * nu); // plane stress
	CC[0][0] = c * 1.;
	CC[0][1] = c * nu;
	CC[0][2] = 0.;
	CC[1][1] = c * 1.;
	CC[1][2] = 0.;
	CC[2][2] = c * 0.5 * (1. - nu);

	for(unsigned int element = 0; element < C.size(); ++element)
		C[element] = utils::tensor::voigt_to_standard_rank4<dim>(CC);
}


inline void make_gaussian_filter(
	std::vector<std::vector<std::pair<size_t, double> > > &gauss_conv_map,
	unsigned int Nx,
	unsigned int Ny,
	double std_blur,
	bool PBC)
{
	std::vector<int> neighbor(2, 0);

	// outside a certain radius, the gaussian weight has a value smaller than min_weight, so we neglect those elements
	double std_blur2 = std_blur * std_blur;
	double min_weight = 0.001;
	assert(min_weight * 2 * 3.14159265 * std_blur2 < 1);
	double radius2 = -std::log(2 * 3.14159265 * std_blur2 * min_weight) * 2. * std_blur2;
	double radius = int(std::sqrt(radius2));

	// get the coordiantes that are in the circle of radius r
	std::set<std::pair<int, int> > gauss_conv_map_coords;

	for(int i = -radius; i <= radius; ++i)
		for(int j = -radius; j <= radius; ++j)
			if(i * i + j * j <= radius2)
				gauss_conv_map_coords.insert(std::make_pair(i, j));

	// use the precalculated relative coordinates to make the neighbours map.
	gauss_conv_map.resize(Nx * Ny);

	for(unsigned int element = 0; element < Nx * Ny; ++element)
	{
		int element_x = element % Nx; // map 1D index into 2D array
		int element_y = std::floor(element / Nx);

		for(auto &p : gauss_conv_map_coords)
		{
			int dx = p.first;
			int dy = p.second;

			neighbor[0] = dx + element_x;
			neighbor[1] = dy + element_y;

			if(PBC)
			{
				neighbor[0] = utils::mod(neighbor[0], Nx);
				neighbor[1] = utils::mod(neighbor[1], Ny);
			}
			else
			{
				if(neighbor[0] < 0 or neighbor[0] >= Nx)
					continue;
				if(neighbor[1] < 0 or neighbor[1] >= Ny)
					continue;
			}

			double weight = std::exp(-double(dx * dx + dy * dy) / (2. * std_blur2));
			unsigned int element_neighbor = neighbor[0] + Nx * neighbor[1];
			gauss_conv_map[element].push_back(std::make_pair(element_neighbor, weight));
		}
	}

	/// ensure that each application of the filter is such that the sum of the weigths is 1
	for(unsigned int element = 0; element < Nx * Ny; ++element)
	{
		double A = 0.0;
		for(auto &p : gauss_conv_map[element])
			A += p.second;

		for(auto &p : gauss_conv_map[element])
			p.second /= A;
	}
}


} // impl
}  // namespace elasticity_solver
} // namespace mepls

#endif //__SolverImpl__