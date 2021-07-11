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

#ifndef RUN_SIM_UTILS_H
#define RUN_SIM_UTILS_H

#include <iostream>
#include <fstream>
#include <random>
#include <math.h>
#include <algorithm>
#include <set>
#include <regex>
#include <string>

#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/timer.h>
#include <deal.II/lac/full_matrix.h>

#if defined(OPENMP)
#include <omp.h>
#endif

#ifndef _COLORS_
#define _COLORS_
	// macros to format strings easily when printed to std::cout
	#define RST  "\x1B[0m"
	#define KRED  "\x1B[31m"
	#define KGRN  "\x1B[32m"
	#define FRED(x) KRED x RST // red
	#define FGRN(x) KGRN x RST // green
	#define BOLD(x) "\x1B[1m" x RST // bold
#endif // _COLORS_


#ifdef DEBUG
#   define M_Assert(Expr, Msg) \
   __M_Assert(#Expr, Expr, __FILE__, __LINE__, Msg)
#else
#   define M_Assert(Expr, Msg) ;
#endif // DEBUG

inline void __M_Assert(const char *expr_str, bool expr, const char *file, int line, const char *msg)
{
	/*! Assert function, which aborts the program execution and prints to
	 * std::cerr a message if the condition is not true. If the preprocessor
	 * variable DEBUG is not defined, the assertion is ignored. */

	// example: M_Assert(ptr != nullptr, "MyFunction: requires non-null argument");

	// from https://stackoverflow.com/a/37264642/5147363

	if(!expr)
	{
		std::cerr << "------------------------------" << "\n" << BOLD(FRED("Assert failed:\t"))
				  << msg << "\n" << "Expected:\t" << expr_str << "\n" << "Source:\t\t" << file
				  << ", line " << line << "\n" << "------------------------------" << "\n";
		abort();
	}
}


namespace mepls
{


/*!
 * This namespce contains tools for testing.
 */
namespace test
{

std::string test_is_passed(bool b)
{
	/*! Return a string formatted for std::cout. If b=true, returns "PASSED" in
	 * bold green; if false, returns "FALSE" in bold red. */

	if(b)
		return BOLD(FGRN("PASSED"));
	else
		return BOLD(FRED("FAILED"));
}

template<int dim>
void assert_result(
	std::string name,
	dealii::SymmetricTensor<4, dim> &actual,
	dealii::SymmetricTensor<4, dim> &expected,
	double tol_percent,
	std::ostream &output_stream = std::cout)
{
	/*! Evaluate if the input "actual" equals the input "expected", up to the
	 * given tolerance (in percentage). A message showing the input name is
	 * printed to output_stream summarizing the result
	 * (see @ref test_is_passed). */

	double z = 100 * (actual - expected).norm() / expected.norm();
	output_stream << " * " << name << " ---- " << z << "% -> " << test_is_passed(z < tol_percent)
				  << std::endl;
}

template<int dim>
void assert_result(
	std::string name,
	dealii::SymmetricTensor<2, dim> &actual,
	dealii::SymmetricTensor<2, dim> &expected,
	double tol_percent,
	std::ostream &output_stream = std::cout)
{
	/*! Evaluate if the input "actual" equals the input "expected", up to the
	 * given tolerance (in percentage). A message showing the input name is
	 * printed to output_stream summarizing the result
	 * (see @ref test_is_passed). */

	double z = 100 * (actual - expected).norm() / expected.norm();
	output_stream << " * " << name << ": val = " << actual << "  |  expected = " << expected << "  |  "
																					"result"
																				 " = "
				  << test_is_passed(z < tol_percent) << std::endl;
}

void assert_result(
	std::string name,
	double actual,
	double expected,
	double tol_percent,
	std::ostream &output_stream = std::cout)
{
	/*! Evaluate if the input "actual" equals the input "expected", up to the
	 * given tolerance (in percentage). A message showing the input name is
	 * printed to output_stream summarizing the result
	 * (see @ref test_is_passed). */

	double z = 100 * std::abs((actual - expected) / expected);
	output_stream << " * " << name << ": val = " << actual << "  |  expected = " << expected <<
	"  |  result"
																				 " = "
				  << test_is_passed(z < tol_percent) << std::endl;
}

} // namespace test


/*!
 * This namespce contains utilities that come in handy in different situations.
 */
namespace utils
{

template<int dim>
inline double get_von_mises_equivalent_stress(const dealii::SymmetricTensor<2, dim> &stress);
/*!< Compute the von Mises equivalent stress of the input stress tensor. */

template<int dim>
inline double get_von_mises_equivalent_strain(const dealii::SymmetricTensor<2, dim> &strain);

/*!< Compute the von Mises equivalent strain of the input strain tensor. */


template<>
inline double get_von_mises_equivalent_stress(const dealii::SymmetricTensor<2, 2> &stress)
{
	/*! Specialization for 2D of the template function @ref
	 * get_von_mises_equivalent_stress(const dealii::SymmetricTensor<2,dim> &)
	 *
	 * @note the von Mises strain \f$ \Sigma_{\rm VM} \f$ is defined here as
	 * \f$ \Sigma_{\rm VM} = \sqrt{\frac{1}{2} \boldsymbol{\Sigma}^{\prime}:
	 * \boldsymbol{\Sigma}^{\prime}} \f$, where \f$ \boldsymbol{\Sigma}^{\prime}
	 * \f$ is the deviatoric part of the stress tensor. Note that alternative
	 * definitions condider an extra factor of \f$ \sqrt{3} \f$. */

	return std::sqrt(
		(1 / 2.) * (0.5 * std::pow(stress[0][0] - stress[1][1], 2.) + 2. * std::pow(stress[0][1],
																					2.)));
}

template<>
inline double get_von_mises_equivalent_strain(const dealii::SymmetricTensor<2, 2> &strain)
{
	/*! Specialization for 2D of the template function @ref
	 * get_von_mises_equivalent_strain(const dealii::SymmetricTensor<2,dim> &)
	 *
	 * @note the von Mises strain \f$ \epsilon_{\rm VM} \f$ is defined here as
	 * \f$ \epsilon_{\rm VM} = \sqrt{2 \boldsymbol{\epsilon}^{\prime}:
	 * \boldsymbol{\epsilon}^{\prime}} \f$, where \f$
	 * \boldsymbol{\epsilon}^{\prime} \f$ is the deviatoric part of the strain
	 * tensor. Note that alternative definitions condider an extra factor of
	 * \f$ 1/\sqrt{3} \f$. */

	return std::sqrt(
		(2.) * (0.5 * std::pow(strain[0][0] - strain[1][1], 2.) + 2. * std::pow(strain[0][1], 2.)));
}



/*! This namespace contains utility functions to work with tensors, mainly to 
 * create and manipulate 4-rank stiffness tensors. */
namespace tensor
{

// the following template functions have a template argument of value 3*(dim-1).
// This is the way to obtain the right dimension for the matrices when using 
// Voigt's notation. Thus, when dim=2, Voigt's matrices are 3x3 and when dim=3, 
// Voigt's matrices are 6x6

template<int dim>
inline dealii::SymmetricTensor<4, dim> voigt_to_standard_rank4(
	dealii::SymmetricTensor<2, 3 * (dim - 1)> const &C);

/*!< Convert a symmetric rank-2 tensor written in Voigt's notation to its
 * symmetric rank-4 form. */

template<>
inline dealii::SymmetricTensor<4, 2> voigt_to_standard_rank4(dealii::SymmetricTensor<2, 3> const &C)
{
	/*! Template specialization of @ref
	 * voigt_to_standard_rank4(dealii::SymmetricTensor<2,3*(dim-1)> const &)
	 * for dim=2. */

	dealii::SymmetricTensor<4, 2> CC;
	CC[0][0][0][0] = C[0][0];
	CC[0][0][1][1] = C[0][1];
	CC[0][0][0][1] = C[0][2];
	CC[1][1][0][0] = C[1][0];
	CC[1][1][1][1] = C[1][1];
	CC[1][1][0][1] = C[1][2];
	CC[0][1][0][0] = C[2][0];
	CC[0][1][1][1] = C[2][1];
	CC[0][1][0][1] = C[2][2];

	return CC;
}

template<int dim>
inline dealii::SymmetricTensor<4, dim> mandel_to_standard_rank4(
	dealii::SymmetricTensor<2, 3 * (dim - 1)> const &C);

/*!< Convert a symmetric rank-2 tensor written in Mandel's notation to its
 * symmetric rank-4 form. */

template<>
inline dealii::SymmetricTensor<4, 2> mandel_to_standard_rank4(
	dealii::SymmetricTensor<2, 3> const &C)
{
	/*! Template specialization of @ref
	 * mandel_to_standard_rank4(dealii::SymmetricTensor<2,3*(dim-1)> const &)
	 * for dim=2. */

	dealii::SymmetricTensor<4, 2> CC;
	CC[0][0][0][0] = C[0][0];
	CC[0][0][1][1] = C[0][1];
	CC[0][0][0][1] = C[0][2] / std::sqrt(2.);
	CC[1][1][0][0] = C[1][0];
	CC[1][1][1][1] = C[1][1];
	CC[1][1][0][1] = C[1][2] / std::sqrt(2.);
	CC[0][1][0][0] = C[2][0] / std::sqrt(2.);
	CC[0][1][1][1] = C[2][1] / std::sqrt(2.);
	CC[0][1][0][1] = C[2][2] / 2.;

	return CC;
}

template<int dim>
inline dealii::SymmetricTensor<2, 3 * (dim - 1)> standard_rank4_to_voigt(
	dealii::SymmetricTensor<4, dim> const &CC);

/*!< Convert a symmetric rank-4 tensor to a symmetric rank-2 one written in
 * Voigt's notation. */

template<>
inline dealii::SymmetricTensor<2, 3> standard_rank4_to_voigt(
	dealii::SymmetricTensor<4, 2> const &CC)
{
	/*! Template specialization of @ref
	 * standard_rank4_to_voigt(dealii::SymmetricTensor<4,2> const &) for dim=2.*/

	dealii::SymmetricTensor<2, 3> C;

	C[0][0] = CC[0][0][0][0];
	C[0][1] = CC[0][0][1][1];
	C[0][2] = CC[0][0][0][1];
	C[1][0] = CC[1][1][0][0];
	C[1][1] = CC[1][1][1][1];
	C[1][2] = CC[1][1][0][1];
	C[2][0] = CC[0][1][0][0];
	C[2][1] = CC[0][1][1][1];
	C[2][2] = CC[0][1][0][1];

	return C;
}


template<int dim>
inline void rotate_matrix(dealii::SymmetricTensor<2, dim> &A, double theta);

/*!< Rotate the input rank-2 tensor an angle of theta radians along the
 * xy-plane. */

template<>
inline void rotate_matrix(dealii::SymmetricTensor<2, 2> &A, double theta)
{
	/*! Template specialization of @ref
	 * rotate_matrix(dealii::SymmetricTensor<2,dim> &, double) for dim=2. */

	dealii::Tensor<2, 2> R;
	R[0][0] = std::cos(theta);
	R[1][1] = R[0][0];
	R[1][0] = std::sin(theta);
	R[0][1] = -R[1][0];

	dealii::Tensor<2, 2> B = R * A * dealii::transpose(R);

	// B is symmetric as well, although dealii returns a normal tensor
	A[0][0] = B[0][0];
	A[1][1] = B[1][1];
	A[0][1] = B[0][1];
}


template<int dim>
inline dealii::SymmetricTensor<2, dim> make_schmid(double theta);

/*!< Create a Schmid tensor (aka slip tensor), defined as
 * \f$ \boldsymbol{M} = \frac{1}{2} \left( \boldsymbol{s} \otimes 
 * \boldsymbol{n} + \boldsymbol{n} \otimes \boldsymbol{s} \right) \f$, where 
 * \f$ \boldsymbol{s} \f$ is the slip direction and \f$ \boldsymbol{n} \f$ is
 * the normal to the slip plane. */


template<>
inline dealii::SymmetricTensor<2, 2> make_schmid<2>(double angle)
{
	/*! Template specialization of @ref make_schmid<dim>(double theta) for dim=2.
	 * Specifically, this function creates a Schmid tensor for a slip plane
	 * forming an angle theta with the horizontal axis. */

	double cos = std::cos(angle);
	double sin = std::sin(angle);

	dealii::SymmetricTensor<2, 2> M;
	M[0][0] = -cos * sin;
	M[1][1] = -M[0][0];
	M[0][1] = 0.5 * (cos * cos - sin * sin);

	return M;
}


template<int dim>
inline dealii::SymmetricTensor<4, dim> make_isotropic_stiffness(double G, double nu);

/*!< Create an isotropic stiffness tensor with shear modulus G and Poisson's
 * ratio nu. */


template<>
inline dealii::SymmetricTensor<4, 2> make_isotropic_stiffness<2>(double G, double nu)
{
	/*! Template specialization of @ref make_isotropic_stiffness(double G,
	 * double nu) for dim=2. */

	dealii::SymmetricTensor<2, 3> CC;

	double E = 2 * G * (1 + nu);
	double c = E / (1. - nu * nu); // plane stress
	CC[0][0] = c * 1.;
	CC[0][1] = c * nu;
	CC[0][2] = 0.;
	CC[1][1] = c * 1.;
	CC[1][2] = 0.;
	CC[2][2] = c * 0.5 * (1. - nu);
	return voigt_to_standard_rank4<2>(CC);
}


template<int dim>
inline dealii::SymmetricTensor<2, 3 * (dim - 1)> make_mandel_anisotropic_stiffness(
	double K, double G1, double G2, double theta);

/*!< Create a stiffness tensor, in Mandel's notation, which is anisotropic to 
 * shear but isotropic to compression.
 * @param K the bulk modulus
 * @param G1 first shear modulus
 * @param G2 second shear modulus
 * @param theta degree of shear anisotropy
 * @note see Nicholas et al., Journal of the Mechanics and Physics of Solids 78 
 * (2015) 333–351 for a detailed explanation. */


template<>
inline dealii::SymmetricTensor<2, 3> make_mandel_anisotropic_stiffness<2>(
	double K, double G1, double G2, double theta)
{
	/*! Template specialization of @ref make_mandel_anisotropic_stiffness(double,
	 * double, double, double) for dim=2. */

	// see in Nicholas et al., Journal of the Mechanics and Physics of Solids 78
	// (2015) 333–351 for further details. Construct a stiffness tensor that is
	// anisotropic to shear but isotropic to volumetric deformation

	// We create it in Mandel notation from the bulk modulus K, and the two shear
	// modules G1 and G2 in the reference frame (rotated by theta), in which it
	// appears as isotropic.

	theta /= 2.;

	// G1 is by definition the smaller of both shear moduli
	double g1 = G1;
	double g2 = G2;
	if(g1 > g2)
	{
		G2 = g1;
		G1 = g2;
	}

	double lambda1 = 2 * G1;
	double lambda2 = 2 * G2;
	double lambda3 = 2 * K;
	double cos2 = std::pow(std::cos(2 * theta), 2.);
	double sin2 = std::pow(std::sin(2 * theta), 2.);
	double a = (lambda3 + lambda1 * sin2 + lambda2 * cos2) / 2.;
	double d = (lambda3 - lambda1 * sin2 - lambda2 * cos2) / 2.;
	double b = (lambda2 - lambda1) * std::sin(4. * theta) * std::sqrt(2.) / 4.;
	double nu = lambda1 * cos2 + lambda2 * sin2;

	dealii::SymmetricTensor<2, 3> C;
	C[0][0] = a;
	C[0][1] = d;
	C[0][2] = b;
	C[1][1] = a;
	C[1][2] = -b;
	C[2][2] = nu;

	return C;
}


template<int dim>
inline dealii::SymmetricTensor<2, 3 * (dim - 1)> compute_voigts_stiffness(
	const dealii::SymmetricTensor<2, dim> &stress_a,
	const dealii::SymmetricTensor<2, dim> &stress_b,
	const dealii::SymmetricTensor<2, dim> &stress_c,
	const dealii::SymmetricTensor<2, dim> &elastic_strain_a,
	const dealii::SymmetricTensor<2, dim> &elastic_strain_b,
	const dealii::SymmetricTensor<2, dim> &elastic_strain_c,
	dealii::Vector<double> &f,
	dealii::Vector<double> &u,
	dealii::FullMatrix<double> &M);

/*!< Compute the stiffness tensor in Voigt's notation from the input elastic
 * strain and stress tensors. The stress and strain tensors denoted as a, b, 
 * and c correspond to 3 different and independent deformation modes (e.g., 
 * tensile xx, tensile yy, shear xy). */


template<>
inline dealii::SymmetricTensor<2, 3> compute_voigts_stiffness(
	const dealii::SymmetricTensor<2, 2> &stress_a,
	const dealii::SymmetricTensor<2, 2> &stress_b,
	const dealii::SymmetricTensor<2, 2> &stress_c,
	const dealii::SymmetricTensor<2, 2> &elastic_strain_a,
	const dealii::SymmetricTensor<2, 2> &elastic_strain_b,
	const dealii::SymmetricTensor<2, 2> &elastic_strain_c,
	dealii::Vector<double> &f,
	dealii::Vector<double> &u,
	dealii::FullMatrix<double> &M)
{
	/*! Template specialization of @ref compute_voigts_stiffness<dim>(...) for
	 * dim=2. */

	// Compute the stiffness tensor CC in Voigt's notation from the input elastic
	// strain and stress in standard notation. The stress and strain a,b,c
	// correspond to 3 different and independent loading modes (e.g., tensile
	// xx, tensile yy, shear xy)

	// we want to compute the components of the rank-4 tensor C. We write the
	// relation between stress (s) and eigenstrain (e) using Voigt's notation as:
	//
	//  |s00|   |C00 C01 C02|   |e00 |
	//  |s11| = |C10 C11 C12| . |e01 |
	//  |s01|  |C20 C21 C22|   |2e02|
	//
	// We have 9 unkowns C_ij (here we don't assume the symmetry C_ij = C_ji,
	// which we impose by symmetrizing the result later). We solve 3 elastic
	// equilibrium problems (a,b,c), each of which with different independent
	// local eigenstrain. From that, we obtain the local stress variation induced
	// by each eigenstrain tensor. We have thus 3 x 3 = 9 stress values and 9
	// eigensrain values. The components of C can be obtained them from the
	// linear system:
	//
	//   |s^a_00|   |e^a_00 e^a_11 e^a_01  0      0      0       0      0      0   |   |C00|
	//   |s^a_11|   |   0     0      0   e^a_00 e^a_11 e^a_01    0      0      0   |   |C01|
	//   |s^a_01|   |   0     0      0     0      0      0    e^a_00 e^a_11 e^a_01 |   |C02|
	//   |s^b_00|   |e^a_00 e^a_11 e^a_01  0      0      0       0      0      0   |   |C10|
	//   |s^b_11|   |   0     0      0   e^a_00 e^a_11 e^a_01    0      0      0   | . |C11|
	//   |s^c_01|   |   0     0      0     0      0      0    e^a_00 e^a_11 e^a_01 |   |C12|
	//   |s^c_00|   |e^a_00 e^a_11 e^a_01  0      0      0       0      0      0   |   |C20|
	//   |s^c_11|   |   0     0      0   e^a_00 e^a_11 e^a_01    0      0      0   |   |C21|
	//   |s^c_01|   |   0     0      0     0      0      0    e^a_00 e^a_11 e^a_01 |   |C22|
	//
	// Once we obtain all the C_ij, we symmetrize it as C_ij <- (C_ij + C_ji) / 2 for i != j.
	// Afterwards, we change C from the rank-2 Voigt's notation to a rank-4 tensor


	f.reinit(9);
	u.reinit(9);
	M.reinit(9, 9);

	f[0] = stress_a[0][0];
	f[1] = stress_a[1][1];
	f[2] = stress_a[0][1];
	f[3] = stress_b[0][0];
	f[4] = stress_b[1][1];
	f[5] = stress_b[0][1];
	f[6] = stress_c[0][0];
	f[7] = stress_c[1][1];
	f[8] = stress_c[0][1];

	M[0][0] = elastic_strain_a[0][0];
	M[0][1] = elastic_strain_a[1][1];
	M[0][2] = elastic_strain_a[0][1] * 2;
	M[1][3] = elastic_strain_a[0][0];
	M[1][4] = elastic_strain_a[1][1];
	M[1][5] = elastic_strain_a[0][1] * 2;
	M[2][6] = elastic_strain_a[0][0];
	M[2][7] = elastic_strain_a[1][1];
	M[2][8] = elastic_strain_a[0][1] * 2;
	M[3][0] = elastic_strain_b[0][0];
	M[3][1] = elastic_strain_b[1][1];
	M[3][2] = elastic_strain_b[0][1] * 2;
	M[4][3] = elastic_strain_b[0][0];
	M[4][4] = elastic_strain_b[1][1];
	M[4][5] = elastic_strain_b[0][1] * 2;
	M[5][6] = elastic_strain_b[0][0];
	M[5][7] = elastic_strain_b[1][1];
	M[5][8] = elastic_strain_b[0][1] * 2;
	M[6][0] = elastic_strain_c[0][0];
	M[6][1] = elastic_strain_c[1][1];
	M[6][2] = elastic_strain_c[0][1] * 2;
	M[7][3] = elastic_strain_c[0][0];
	M[7][4] = elastic_strain_c[1][1];
	M[7][5] = elastic_strain_c[0][1] * 2;
	M[8][6] = elastic_strain_c[0][0];
	M[8][7] = elastic_strain_c[1][1];
	M[8][8] = elastic_strain_c[0][1] * 2;

	// solve linear system
	M.invert(M);
	M.vmult(u, f);

	// map the vector of unkowns to its tensor form (in Voigt's notation)
	dealii::Tensor<2, 3> CC;
	CC[0][0] = u[0];
	CC[0][1] = u[1];
	CC[0][2] = u[2];
	CC[1][0] = u[3];
	CC[1][1] = u[4];
	CC[1][2] = u[5];
	CC[2][0] = u[6];
	CC[2][1] = u[7];
	CC[2][2] = u[8];

	// symmetrize the result
	dealii::SymmetricTensor<2, 3> symC;
	symC[0][1] = (CC[0][1] + CC[1][0]) / 2.;
	symC[0][2] = (CC[0][2] + CC[2][0]) / 2.;
	symC[1][2] = (CC[1][2] + CC[2][1]) / 2.;
	symC[0][0] = CC[0][0];
	symC[1][1] = CC[0][0];
	symC[2][2] = CC[2][2];

	return symC;
}

} // namespace tensor



/*! This namespace contains utility functions to generate random numbers 
 * according to specific probability distributions. */
namespace rand
{

inline std::vector<std::vector<double>> Cholesky_decomposition(std::vector<std::vector<double>> &matrix)
{
	/*! Compute the Cholesky decomposition of the input matrix. */

	// adapted from https://www.geeksforgeeks.org/cholesky-decomposition-matrix-decomposition/

	unsigned int n = matrix.size();

	std::vector<std::vector<double>> lower(n, std::vector<double>(n, 0.));

	// Decomposing a matrix into Lower Triangular
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j <= i; j++)
		{
			double sum = 0;

			if(j == i) // summation for diagnols
			{
				for(int k = 0; k < j; k++)
					sum += pow(lower[j][k], 2);
				lower[j][j] = std::sqrt(matrix[j][j] - sum);
			}
			else
			{

				// Evaluating L(i, j) using L(j, j)
				for(int k = 0; k < j; k++)
					sum += (lower[i][k] * lower[j][k]);
				lower[i][j] = (matrix[i][j] - sum) / lower[j][j];
			}
		}
	}

	return lower;
}


template<int dim>
inline std::vector<std::vector<double>> draw_multivariate_gaussian_sample(
	std::vector<double> &mu,
	dealii::SymmetricTensor<2, 3 * (dim - 1)> &sym_cov,
	unsigned int N,
	std::mt19937 &generator);

/*!< Draw a sample of size N of multivariate gaussian data points. The aim of
 * this function is to generate local stress or strain tensors, which are 
 * symmetric. Thus, the dimension of the points is 3*(dim-1) since for dim=2,
 * there are three different symmetric tensor components, while for dim=3,
 * there are 6.
 * 
 * @param mu average
 * @param sym_cov covariance matrix
 * @param N size of the sample
 * @returns a vector of 3*(dim-1) vectors, each one of size N. */

template<>
inline std::vector<std::vector<double>> draw_multivariate_gaussian_sample<2>(
	std::vector<double> &mu,
	dealii::SymmetricTensor<2, 3> &sym_cov,
	unsigned int N,
	std::mt19937 &generator)
{
	/*! Template specialization of @ref
	 * draw_multivariate_gaussian_sample<dim>(...) for dim=2. */

	int d = mu.size();

	assert(d == 3);

	std::vector<std::vector<double>> cov(3, std::vector<double>(3));
	cov[0][0] = sym_cov[0][0];
	cov[1][1] = sym_cov[1][1];
	cov[2][2] = sym_cov[2][2];
	cov[0][1] = sym_cov[0][1];
	cov[1][0] = sym_cov[1][0];
	cov[0][2] = sym_cov[0][2];
	cov[2][0] = sym_cov[2][0];
	cov[1][2] = sym_cov[1][2];
	cov[2][1] = sym_cov[2][1];

	// compute and use the Cholesky decomponsition of the covariance matrix to
	// generate multivariate gaussian samples from a univariate distribution
	std::vector<std::vector<double>> C = Cholesky_decomposition(cov);
	std::normal_distribution<double> normal_dist;

	std::vector<std::vector<double>> Y;
	std::vector<double> y(d, 0.);
	std::vector<double> x(d, 0.);
	for(unsigned int n = 0; n < N; ++n)
	{
		for(unsigned int i = 0; i < d; ++i)
			x[i] = normal_dist(generator);

		for(unsigned int i = 0; i < d; ++i)
		{
			y[i] = mu[i];

			for(unsigned int j = 0; j < d; ++j)
				y[i] += C[i][j] * x[j];
		}

		Y.push_back(y);
	}

	return Y;
}

inline double get_exponential_rand(double av, double rand)
{
	/*! Generate an exponentially distributed random number.
	 * @param av average
	 * @param rand uniformly random number distributed between 0 and 1. */

	return -std::log(rand) * av;
}

inline double get_weibull_rand(double k, double l, double rand)
{
	/*! Generate a Weibull distributed random number.
	* @param k shape parameter
	* @param l scale parameter
	* @param rand uniformly random number distributed between 0 and 1. */

	return l * std::pow(-std::log(rand), 1. / k);
}

inline double get_weibull_rand_av(double k, double av, double rand)
{
	/*! Generate a Weibull distributed random number.
	* @param k shape parameter
	* @param av average
	* @param rand uniformly random number distributed between 0 and 1. */

	double l = av / std::tgamma(1. + 1. / k);
	return l * std::pow(-std::log(rand), 1. / k);
}


template<typename FuncType>
std::piecewise_linear_distribution<double> *create_distribution(
	FuncType &function, double dx, double x0, double wmin, double x1)
{
	/*! Create a probability distribution using the input function. It creates a
	 * support interval starting at x0, and in steps of dx the support interval
	 * is increased until reaching a maximum value x1 or until the probability
	 * at x is smaller than wmin.
	 *
	 * @param function the functional form of the pdf
	 * @param dx discretization step of the random variable
	 * @param x0 mimnimum value of x
	 * @param x1 maximum value of x
	 * @param wmin tolerance for deciding wheter the probability of [x,x+dx) is 0
	 *
	 * @return a pointer to a distribution that can be called as
	 * distribution(generator) to generate a random number with a pdf given by
	 * the input function.
	 *
	 * @warning the user is reponsible for deleting the dynamically-allocated
	 * distribution.
	 *
	 * Example of use:
	   @code
		   // define the funcitonal form of the pdf
		  auto weibull = [&](double x, double k, double lambda) { return std::pow(x,k-1)*std::exp(-std::pow(x/lambda,k)); };
		  // set the parameters to specific values, in this case k=2.3 and labmda=4.2. The only free parameter is the random variable x
		  auto func = std::bind(weibull, std::placeholders::_1, 2.3, 4.2);
		  // get a pointer to the callable pdf
		  auto threshold_distribution_ptr = utils::create_random_dist(func, 1e-3, 0., 1e-7, 8);
	   @endcode
	*/

	std::vector<double> intervals;
	std::vector<double> weights;
	double x = x0;
	double w = 1.;

	while(w > wmin or x < x1)
	{
		x += dx;
		intervals.push_back(x);
		w = function(x);
		weights.push_back(w);
	}

	auto dist = new std::piecewise_linear_distribution<double>(intervals.begin(), intervals.end(),
															   weights.begin());

	return dist;
}


} // namespace rand



/*! This namespace contains utility functions to parse and generate strings. */
namespace str
{

inline std::string to_string(const bool b)
{
	/*! Return a string with the value of the input bool, i.e., "true" if b==true
	 * and "false" if b==0. */

	std::ostringstream ss;
	ss << std::boolalpha << b;
	return ss.str();
}

template<typename T>
inline std::string to_string(const T x, const int n = 2)
{
	/*! Return a string with the value of the input parameter of type T. The
	 * string contains the value up to n decimal digits. */

	std::ostringstream out;
	out.precision(n);
	out << std::fixed << x;
	return out.str();
}

std::vector<int> parse_list_integers(const std::string &str)
{
	/*! Parse a string containing a comma-separated list of integers. The values
	 * are returned as a vector. */

	// from https://stackoverflow.com/a/43209178/5147363

	typedef std::regex_iterator<std::string::const_iterator> re_iterator;
	typedef re_iterator::value_type re_iterated;

	std::regex re("(\\d+)");

	re_iterator rit(str.begin(), str.end(), re);
	re_iterator rend;

	std::vector<int> result;

	std::transform(rit, rend, std::back_inserter(result), [](const re_iterated &it)
	{ return std::stoi(it[1]); });

	return result;
}


void read_csv(std::vector<double> &xData, std::vector<double> &yData, std::string filename,
				unsigned int skip=0)
{
	/*! Read a text file that contains two columns of numerical values. The columns
	 * are parsed into two vectors of doubles.
	 *
	 * @param xData vector to store column #1
	 * @param yData vector to store column #2
	 * @param filename path tot the file
	 * @param skip number of rows to skip
	 * */

	// adapted from https://stackoverflow.com/a/30181684

	std::ifstream ifile(filename);

	std::string line; // we read the full line here

	// skip these first lines
	for(unsigned int n=0; n < skip; ++n)
		std::getline(ifile, line);

	while (std::getline(ifile, line)) // read the current line
	{
		std::istringstream iss{line}; // construct a string stream from line

		// read the tokens from current line separated by comma
		std::vector<std::string> tokens; // here we store the tokens
		std::string token; // current token
		while (std::getline(iss, token, ','))
			tokens.push_back(token); // add the token to the vector

		xData.push_back( std::stof(tokens[0]) );
		yData.push_back( std::stof(tokens[1]) );
	}

}

} // namespace string



/*! This class controls whether a simulation should continue running or should 
 * stop. When the object is called using the call operator with some input 
 * condition, if the condition evaluates to false, the object changes its 
 * internal state to false. That state is not changed anymore during the 
 * object's life. When the object is called without input conditions, we query
 * the object's internal state for whether the simulation should continue or 
 * note. The advantage is that we can use the object to evaluate many stopping
 * criteria at different parts of the code, and if one of them is not met,
 * the simulation will stop, and the object will report with an informative 
 * message, which was the specific unfulfilled condition. */
class ContinueSimulation
{
  public:
	ContinueSimulation()
	{
		/*! Constructor. */

		continue_bool = true;
		message = " continue_sim print";
	};

	bool operator()()
	{
		/*! The returned bool tells whether the simulation must continue or
		 * stop. */

		return continue_bool;
	};

	void operator()(bool b, std::string const input_message)
	{
		/*! Set the state of the object to b. If the state was already false, do
		 * nothing. If the state is set to false in this call, the input message
		 * is stored so it can be used later to inform what was the reason for
		 * setting the state to false. */

		// once continue==false, it can't be changed to true
		if(not continue_bool)
			return;

		continue_bool = b;
		message = input_message;
	};

	friend std::ostream &operator<<(std::ostream &stream, const ContinueSimulation &continue_sim)
	{
		/*! Defining the stream insertion operator allows us to print the message
		 * that was given to the object when its state was set to false. */

		stream << "---- " << continue_sim.message << " ----" << std::endl;

		return stream;
	};

  private:
	bool continue_bool;
	/*!< Bool that defines the internal state of the object. If it is true, the
	 * simulation must continue running. If it is false, the simulation must
	 * stop. */

	std::string message;
	/*!< Message associated with the condition that failed. The failed condition
	 * is the reason why the simulation must stop. */
};


template<typename T>
void remove_duplicates(std::vector<T> &vec)
{
	/*! Remove duplicate values in the input vector (i.e., keep only unique
	 * values). */

	// from https://stackoverflow.com/a/9237226/5147363

	std::sort(vec.begin(), vec.end());
	vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
}


inline double mod(double a, double b)
{
	/*! Compute the value of a with modulus b. Example: mod(3,5) = 3,
	 * mod(-1,5) = 4. This definition is useful for performing operation with
	 * periodic boundaries. */

	// from https://stackoverflow.com/a/12277233/228965]

	return std::fmod(std::fmod(a, b) + b, b);
}



/*! Singleton class to measure execution times. The advantage of using a
 * singleton pattern is that we can measure the execution times of any part of
 * the code without the need to modify function interfaces. Specifically, in
 * every call to a function, we can get a pointer to an existing instance of
 * the timer object, avoiding creating a new one. Moreover, since the object
 * is always the same one, the execution times add up, which provides us with
 * an averaged execution time should the call to the function be repeated 
 * iteratively. */
class TimerSingleton
{
  public:

	static TimerSingleton *getInstance()
	{
		/*! Get a pointer to the singleton timer. */

		if(!instanceFlag)
		{
			single = new TimerSingleton();
			instanceFlag = true;
			return single;
		}
		else
		{
			return single;
		}
	}

	void enter_subsection(std::string name)
	{
		/*! Enter a measured subsection with the input name. */

		timer->enter_subsection(name);
	}

	void leave_subsection(std::string name)
	{
		/*! Leave the measured subsection with the input name. */

		timer->leave_subsection(name);
	}

	void print_summary()
	{
		/*! Print the summary of execution times. */

		timer->print_summary();
	}

	~TimerSingleton()
	{
		/*! Desctuctor. */

		instanceFlag = false;
		delete timer;
	}


  private:
	static bool instanceFlag;
	/*!< Bool indicating whether the singleton object has already been created
	 * or not. */

	static TimerSingleton *single;
	/*!< Pointer to the singleton object. */

	dealii::TimerOutput *timer;

	/*!< Pointer to a dealII timer object, providing the timer functionality to
	 * the class. */

	TimerSingleton()
	{
		/*! Private constructor. */

		timer = new dealii::TimerOutput(std::cout, dealii::TimerOutput::summary,
										dealii::TimerOutput::wall_times);
	}

};

bool TimerSingleton::instanceFlag = false;

TimerSingleton *TimerSingleton::single = NULL;



/*! This class provides linear interpolations of a 1D function, defined by the
 * given x and y data vectors. */
class Interpolator
{
  public:

	Interpolator(std::vector<double> &xData_, std::vector<double> &yData_)
		:
	xData(xData_),
	yData(yData_)
	{
		/*! Constructor. */
	}

	double operator()(double x)
	{
		/*! Returns the value of the function, interpolated at the input value x. */

		// adapted from http://www.cplusplus.com/forum/general/216928/

	   int size = xData.size();

		// find left end of interval for interpolation
	   int i = 0;
	   // special case: beyond right end
	   if ( x >= xData[size - 2] )
	   {
		  i = size - 2;
	   }
	   else
	   {
		  while ( x > xData[i+1] ) i++;
	   }
	   // points on either side (unless beyond ends)
	   double xL = xData[i], yL = yData[i], xR = xData[i+1], yR = yData[i+1];

 	   if ( x < xL ) yR = yL;
  	   if ( x > xR ) yL = yR;

		// gradient
	   double dydx = ( yR - yL ) / ( xR - xL );

 		// linear interpolation
	   return yL + dydx * ( x - xL );
	}

  private:

	std::vector<double> xData;
	/*!< Discrete x values of the interpolated function. */

	std::vector<double> yData;
	/*!< Discrete y values of the interpolated function. */
};


} // namespace utils

} // namespace mepls

#endif //RUN_SIM_UTILS_H
