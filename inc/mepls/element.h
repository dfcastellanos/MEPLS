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

#ifndef RUN_SIM_ELEMENTS_H
#define RUN_SIM_ELEMENTS_H

#include <deal.II/base/symmetric_tensor.h>

#include <mepls/utils.h>

namespace mepls
{

namespace element
{
template<int dim>
class Element;

template<int dim> using Vector = typename std::vector<Element<dim> *>;
}

namespace event
{
template<int dim>
class RenewSlip;

template<int dim>
class Plastic;
}

namespace elasticity_solver
{
template<int dim>
class Solver;
}

/*! This namespace contains the abstract slip class slip::Slip, which servers as
 * the interface for slip objects. Slip objects contain the physics of the
 * activation of slip system within mesoscale elements (see \ref element). */
namespace slip
{

/*! Slip objects contain the physics of the activation of the slip system within
 * mesoscale elements. This abstract class provides an interface for the slip
 * objects and specific activation rules are implemented in derived classes.
 * Essentially, slip systems are defined by a slip plane, an effective shear
 * stress on that plane, and the critical effective shear stress known as the
 * threshold that defines when the system becomes mechanically activated. The
 * interplay between these quantities is implemented by derived classes.
 *
 * @note this class defines only the slip activation rules but **does not**
 * represent an active slip event, which is represented by event::Plastic.
 */
template<int dim>
class Slip
{
  public:

	Slip()
		:
		eff_shear_stress(0.),
		pressure(0.),
		threshold(0.),
		barrier(0.),
		activation_rate(0.),
		angle(0.)
	{
		/*! Constructor. Initialize members zero. */
	}

	virtual ~Slip()
	{
		/*! Virtual destructor. */
	}

	virtual void update() = 0;
	/*!< Update the state of the slip object using on the state of the parent
	 * mesoscale element, \ref parent. */

	virtual double get_critical_load_increment() = 0;
	/*!< Calculate the necessary external load increment needed to make the slip
	 * activation barrier vanish. */

	virtual dealii::SymmetricTensor<2, dim> get_eigenstrain_increment() = 0;

	/*!< Calculate the plastic strain increment tensor that represents, at the
	 * mesoscale, the local deformation induced by the slip event. */

	virtual Slip<dim> *make_copy()
	{
		/*! Return a pointer to a dynamically allocated copy of the slip object.
		 * The dynamic type of the object is that of the derived class.
		 *
		 * @warning the user is responsible for deleting the allocated copy. */

		M_Assert(false, "make_copy() not implemented for the current slip object.");
	}

	double eff_shear_stress;
	/*!< Effective shear stress \f$ \tau \f$ in the slip system. It is
	 * responsible for lowering the slip activation barrier, \ref barrier. */

	double pressure;
	/*!< Hydrostatic pressure, computed from the parent element stress tensor
	 * element::\ref stress_ as
	 * \f$ p = \frac{1}{dim}\textrm{tr}( \boldsymbol{\Sigma} )\f$.
	 *
	 * @note the pressure is the same for all the slip objects contained within
	 * the same parent element, since the pressure is anisotropic and only a
	 * single stress tensor per element is considered. */

	double threshold;
	/*!< Local yield threshold \f$ \hat{\tau} \f$, i.e. the effective shear
	 * stress \ref eff_shear_stress necessary for a slip system to become
	 * mechanically unstable. */

	double barrier;
	/*!< Stress barrier \f$ \hat{\tau} - \tau \f$ neccesary to be overcome
	 * before the slip system becomes mechanically activated. */

	double activation_rate;
	/*!< Rate \f$ \nu \f$ for thermal slip activation. */

	double angle;
	/*!< Orientation of the slip plane (in radians). */

	element::Element<dim> *parent;
	/*!< Parent mesoscale element element::Element to which the slip object
	 * belongs. When the state of the slip oject is updated via \ref update(),
	 * this element is accessed.
	 *
	 * @warning parent elements contain pointers to the owned slip objects, thus
	 * deletion of slips objects must be done carefully. Usually, allocation and
	 * deletion of slip objects is only done within the implementation of
	 * element::Element objects. */
};

} // namespace slip


/*! Element objects represent mesoscale regions of a material. This namespace
 * contains the abstract element class element::Element, which serves as the
 * interface for different element objects. It also contains free functions
 * extending its interface */
namespace element
{

/*! This struct carries the information about how the microstructural properties
 * must be renewed within a mesoscale element after a slip event has taken
 * place within that element. */
template<int dim>
struct RenewInstruct
{
	bool slip_properties = true;
	/*!< Renew the slip systems owned by the element. */

	bool elastic_properties = false;
	/*!< Renew the elastic properties of the element. */

	event::Plastic<dim> plastic_event;
	/*!< Plastic event responsible of the change of microstructural properties.
	 * Its details, such as, e.g., the amplitude of the plastic deformation, can
	 * be taken into account during the renewal process. */
};


/*! Element objects represent local mesoscale regions of a material. The state
 * of an element is defined by continuum mechanics magnitudes (stress, strain
 * and eigenstrain tensors) and slip systems present in the represented mesoscale
 * region. Each element has associated a single tensor for each magnitude. The
 * behavior of the slip systems (defined by objects of class slip::Slip) define
 * the plastic response of an element. The local element microstructural
 * properties are defined by the element's elastic properties and the existing
 * slip systems. Although from a physical perspective, each mesocale element has
 * a well-defined spatial location, the element object here defined is not aware
 * about spatial coordinates. This is so because an element's behavior is
 * defined entirely in terms of local fields and slip properties. Thus, a
 * material discretized in a lattice of mesoscale elements can be represented
 * as a 1D vector of element objects, while spatial information is present only
 * during the  computation of the elastic fields, which is performed outside
 * the element objects.
 *
 * This abstract class provides an interface for element objects and specific
 * implementations are defined by derived classes. */

template<int dim>
class Element
{
  public:

	Element()
	{
		/*! Constructor. Initialize the default values of the class members. */
		integrated_vm_eigenstrain_ = 0.;
		S_already_is_set = false;
		ext_stress_coeff_is_set = false;
		this->type("element");
	}

	Element(Element *element)
	{
		this->C( element->C() );
		this->S( element->S());
		this->eigenstrain( element->eigenstrain() );
		this->integrated_vm_eigenstrain( element->integrated_vm_eigenstrain() );
		this->ext_stress_coeff( element->ext_stress_coeff() );
		this->prestress( element->prestress() );
		this->elastic_stress( element->elastic_stress() );
		this->def_grad( element->def_grad() );
		this->type( element->type() );
		this->number( element->number() );

		for(auto &slip : *element)
			this->add_slip_system( slip->make_copy() );
	}

	virtual ~Element()
	{
		/*! Virtual destructor */

		// Call remove_slip_systems to remove the dynamically allocated slip
		// objects owned by the element.
		remove_slip_systems();
	}

	virtual void record_structural_properties(std::vector<event::RenewSlip<dim>> &renewal_vector)
	{
		/*! Write the current structural properties of the element into the
		 * input object.
		 *
		 * @note at the moment, only slip properties are written. */

		for(auto &slip : *this)
		{
			event::RenewSlip<dim> renewal_event(slip);
			renewal_vector.push_back(renewal_event);
		}
	}

	void renew_structural_properties(RenewInstruct<dim> &renew_instruct)
	{
		/*! Renew the structural properties of the element (e.g., renewing slip
		 * thresholds, slip angles, elastic properties, etc.). Specific renewal
		 * rules are defined in derived classes by implementing \ref
		 * renew_structural_properties_impl. The input RenewInstruct object
		 * carries information that might be taken into account when renewing
		 * the properties. */

		renew_structural_properties_impl(renew_instruct);
	}

	void renew_structural_properties()
	{
		/*! Overloaded version of \ref renew_structural_properties which, for
		 * conviniency, takes no input. It wraps the original function and
		 * passes to it a default-initialised RenewInstruct object. */

		RenewInstruct<dim> renew_instruct;
		renew_structural_properties(renew_instruct);
	}

	virtual void renew_structural_properties_impl(RenewInstruct<dim> &renew_instruct) = 0;

	/*!< Define how the element slip and elastic properties are to be renewed
	 * when \ref renew_structural_properties is called. This abstract function
	 * is to be implemented by derived classes. */

	mepls::element::Element<dim> *make_copy()
	{
		/*! Return a pointer to a dynamically allocated copy of the object. The
		 * dynamic type of the object is that of the derived class.
		 *
		 * @warning the user is responsible for deleting the allocated copy. */

		return make_copy_impl();
	}

	void add_eigenstrain(const dealii::SymmetricTensor<2, dim> &eigenstrain_increment)
	{
		/*! Add an eigenstrain increment
		 * \f$ \Delta\boldsymbol{\epsilon}_{\rm pl} \f$ to the element. Thus,
		 * the element's eigenstrain \ref eigenstrain_ is updated as
		 * \f$ \boldsymbol{\epsilon}_{\rm pl} +=
		 * \Delta \boldsymbol{\epsilon}_{\rm pl}\f$ and the integrated von Mises
		 * eigenstrain \ref integrated_vm_eigenstrain_ as \f$ \epsilon_{\rm VM}
		 * += \sqrt{\frac{1}{2} \Delta \boldsymbol{\epsilon}^{\prime}_{\rm pl}:
		 * \Delta \boldsymbol{\epsilon}^{\prime}_{\rm pl}} \f$, where primed
		 * magnitudes refer to deviatoric tensors. */

		eigenstrain_ += eigenstrain_increment;
		integrated_vm_eigenstrain_ += utils::get_von_mises_equivalent_strain(eigenstrain_increment);
	}

	double integrated_vm_eigenstrain() const
	{
		/*! Return the value of the integrated von Mises eigenstrain
		 * \ref integrated_vm_eigenstrain_. */

		return integrated_vm_eigenstrain_;
	}

	const dealii::SymmetricTensor<2, dim> &eigenstrain() const
	{
		/*! Return a reference to the accumlated eigenstrain tensor
		 * \ref eigenstrain_. */

		return eigenstrain_;
	}

	const dealii::SymmetricTensor<2, dim> &prestress() const
	{
		/*! Return a reference to the prestress tensor \ref presstress_. */

		return prestress_;
	}

	void prestress(const dealii::SymmetricTensor<2, dim> &prestress)
	{
		/*! Set a new prestress tensor \ref presstress_. The total stress is
		 * updated afterwards. */

		prestress_ = prestress;

		update_stress();
	}

	const dealii::SymmetricTensor<2, dim> &stress() const
	{
		/*! Return a reference to the total stress tensor \ref stress_. */

		return stress_;
	}

	const dealii::SymmetricTensor<2, dim> &elastic_stress() const
	{
		/*! Return a reference to the elastic stress tensor
		 * \ref elastic_stress_. */

		return elastic_stress_;
	}

	void elastic_stress(const dealii::SymmetricTensor<2, dim> &input_elastic_stress)
	{
		/*! Set a new elastic stress tensor \ref elastic_stress_. The total
		 * stress is updated afterwards as defined
		 * by \ref update_stress. */

		elastic_stress_ = input_elastic_stress;

		update_stress();
	}

	const dealii::Tensor<2, dim> &def_grad() const
	{
		/*! Return a reference to the deformation grandient tensor,
		 * \ref def_grad_. */

		return def_grad_;
	}

	void def_grad(const dealii::Tensor<2, dim> &input_def_grad)
	{
		/*! Set a new deformation gradient tensor, \ref def_grad_. */

		def_grad_ = input_def_grad;
	}

	void S(const dealii::SymmetricTensor<4, dim> &input_S)
	{
		/*! Set a new \ref S_ tensor. */

		S_ = input_S;
		S_already_is_set = true;
	}

	const dealii::SymmetricTensor<4, dim> &S() const
	{
		/*! Return a reference to the \ref S_ tensor. */

		M_Assert(S_already_is_set, "Local strain coeffs. have not been set yet");

		return S_;
	}

	void C(const dealii::SymmetricTensor<4, dim> &input_C)
	{
		/*! Set a new \ref C_ tensor. */

		C_ = input_C;

		update_stress();
	}

	const dealii::SymmetricTensor<4, dim> &C() const
	{
		/*! Return a reference to the \ref C_ tensor. */

		return C_;
	}

	const dealii::SymmetricTensor<2, dim> &ext_stress_coeff() const
	{
		/*! Return a reference to the \ref ext_stress_coeff_ tensor. */

		M_Assert(ext_stress_coeff_is_set, "Ext. stress coeffs. have not been set yet");

		return ext_stress_coeff_;
	}

	void ext_stress_coeff(const dealii::SymmetricTensor<2, dim> &ext_stress_coeff)
	{
		/*! Set a new \ref ext_stress_coeff_ tensor. */

		ext_stress_coeff_ = ext_stress_coeff;
		ext_stress_coeff_is_set = true;
	}

	typename std::vector<mepls::slip::Slip<dim> *>::iterator begin()
	{
		/*! Beginning of the iteration over the element. It corresponds to the
		 * beginning of the vector containing the slip objects owned by the
		 * element. */

		return slip_systems_.begin();
	}

	typename std::vector<mepls::slip::Slip<dim> *>::iterator end()
	{
		/*! End of the iteration over the element. It corresponds to the end of
		 * the vector containing the slip objects owned by the element. */

		return slip_systems_.end();
	}

	void add_slip_system(slip::Slip<dim> *slip)
	{
		/*! Add a new slip object to the element.
		 *
		 * @warning the slip objects added to an element must be dynamically
		 * allocated, since they will be deleted when calling
		 * \ref remove_slip_systems. */

		// inform the slip object about its parent
		slip->parent = this;

		// update the slip object state so it takes into account its parent's
		// stress state
		slip->update();

		// add the slip object to the element
		slip_systems_.push_back(slip);
	}

	void remove_slip_systems()
	{
		/*! Delete the dynamically allocated slip objects and clear the
		 * vector containing the pointers to the deleted objects. */

		for(slip::Slip<dim> *slip : slip_systems_)
			delete slip;

		slip_systems_.clear();
	}

	slip::Slip<dim> *operator[](unsigned int n)
	{
		/*! Return the slip system number n, as stored in \ref slip_systems_.
		 * It is recommended to accesses slip systems by iterating over the
		 * element instead of indexing. */

		M_Assert(n < slip_systems_.size(), "");

		return slip_systems_[n];
	}

	unsigned int size() const
	{
		/*! Return the size of the element, defined as the number of slip
		 * objects that it currently owns. */

		return slip_systems_.size();
	}

	unsigned int number() const
	{
		/*! Return the number of the element. */

		return number_;
	}

	void number(unsigned int number)
	{
		/*! Set the number of the element. */

		number_ = number;
	}

	void type(const std::string &type)
	{
		/*! Set the type of the element. */

		type_ = type;
	}

	const std::string &type() const
	{
		/*! Return the type of the element. */

		return type_;
	}

	void state_to_prestress()
	{
		/*! Set the current total stress as prestress and clean the rest of
		 * the deformation history, i.e. the eigenstrain, integrated von Mises
		 * eigenstrain, elastic stress, and deformation gradient. */

		prestress_ = stress_;

		eigenstrain_.clear();
		integrated_vm_eigenstrain_ = 0.;
		elastic_stress_.clear();
		stress_.clear();
		def_grad_.clear();

		update_stress();
	}

  private:

	dealii::SymmetricTensor<2, dim> eigenstrain_;
	/*!< Eigenstrain \f$ \boldsymbol{\epsilon}_{\rm pl} \f$ of the element. It
	 * is the acumulation of the eigenstrain increments added with
	 * \ref add_eigenstrain(). */

	double integrated_vm_eigenstrain_;
	/*!< The sum of von Mises eigenstrain increments. The increments are
	 * computed from each eigenstrain tensorial increment added with
	 * \ref add_eigenstrain(). */

	dealii::SymmetricTensor<2, dim> prestress_;
	/*!< Stress \f$ \boldsymbol{\Sigma}_{\rm 0} \f$ that the element has prior
	 * to loading the material. Physically, this field is ypically originated
	 * due to prior plastic deformation or during the solidification of a melt
	 * from where the solid material is cast. */

	dealii::SymmetricTensor<2, dim> elastic_stress_;
	/*!< Elastic stress, \f$ \boldsymbol{\Sigma}_{\rm el} \f$. Its computation
	 * is performed outside the element objects, i.e., elements do not care
	 * about its origing. Therefore, altouth it depends on how the stress fields
	 * are computed, typically this stress is computed from the elastic strain
	 * by linear elasticity, and is the superposition of the internal stress
	 * field consequence of plastic strain, and the external stress field
	 * induced by an external loading mechanism. */

	dealii::SymmetricTensor<2, dim> stress_;
	/*!< Total stress \f$ \boldsymbol{\Sigma} = \boldsymbol{\Sigma}_{\rm 0} +
	 * \boldsymbol{\Sigma}_{\rm el}\f$. It is the superposition of the elastic
	 * stress and the presstress. */

	dealii::Tensor<2, dim> def_grad_;
	/*!< Deformation graient tensor \f$ \boldsymbol{F} \f$. */

	dealii::SymmetricTensor<4, dim> S_;
	/*!< 4-rank tensor which relates local elastic stress variations induced by
	 * local eigesntrain variations. Thus, \f$ \boldsymbol{\Sigma}_{\rm el} =
	 * \mathbb{S} : \Delta \boldsymbol{\epsilon}_{\rm pl} \f$. It is related
	 * to the stiffness \f$ \mathbb{C} \f$ and Eshelby's \f$ \mathbb{E} \f$
	 * 4-rank tensors as \f$ \mathbb{S} = \mathbb{C} :
	 * (\mathbb{E}-\mathbb{I}) \f$, where \f$ \mathbb{I} \f$ is the rank-4
	 * identity tensor. */

	dealii::SymmetricTensor<4, dim> C_;
	/*!< rank-4 stiffness tensor \f$ \mathbb{C} \f$ characterising the element's
	 * linear elastic response, i.e., \f$ \boldsymbol{\Sigma}_{\rm el} =
	 * \mathbb{C} : \boldsymbol{\epsilon}_{\rm el} \f$, where : refers to the
	 * double contraction product. */

	dealii::SymmetricTensor<2, dim> ext_stress_coeff_;
	/*!< Local elastic stress tensor increment \f$ \boldsymbol{a} \f$, induced
	 * by an external load increment of unit amplitude. */

	std::vector<slip::Slip<dim> *> slip_systems_;
	/*!< Vector containing the pointers to the slip objects owned by the
	 * element. */

	unsigned int number_;
	/*!< Element number. */

	std::string type_;
	/*!< String containing the type of element object.*/

	bool S_already_is_set;
	/*!< Indicates wheter the \ref S_ tensor has already been set. */

	bool ext_stress_coeff_is_set;

	/*!< Indicates whether the ext_stress_coeff_ tensor has already been set. */

	void update_stress()
	{
		/*! Compute the total stres, which is a superposition of the prestress
		 * and the elastic stress, i.e., \f$ \boldsymbol{\Sigma} =
		 * \boldsymbol{\Sigma}_{\rm 0} +
		 * \boldsymbol{\Sigma}_{\rm el}\f$. Afterwards, update the stress state
		 * of the slip systems owned by the element. */

		stress_ = prestress_ + elastic_stress_;

		for(auto &slip : *this)
			slip->update();
	}

	void eigenstrain(const dealii::SymmetricTensor<2, dim> &eigenstrain)
	{
		/*! Set a new eigenstrain tensor. Since the user controls the
		 * eigenstrain only through \ref add_eigenstrain, this function is only
		 * used for setting the state of the eigenstrain member when an element
		 * object is copied. */

		eigenstrain_ = eigenstrain;
	}

	void integrated_vm_eigenstrain(double input_integrated_vm_eigenstrain)
	{
		/*! Set a new value for integrated_vm_eigenstrain_. Since the user
		 * controls the eigenstrain only through \ref add_eigenstrain, this
		 * function is only used for setting the state of the
		 * integrated_vm_eigenstrain_ member when a element object is copied. */

		integrated_vm_eigenstrain_ = input_integrated_vm_eigenstrain;
	}

	virtual mepls::element::Element<dim> *make_copy_impl()
	{
		/*! Dynamically allocate a new object of the derive class type and
		 * return a pointer to it. */

		M_Assert(false, "Copying of the current type of element not implemented.");

		abort();
	}
};


template<int dim>
void calculate_ext_stress_coefficients(
	element::Vector<dim> &elements, elasticity_solver::Solver<dim> &solver)
{
	/*! Calculate, using the input elasticity solver, the stress variation
	 * induced in each element by an external load increment of amplitude 1.0.
	 * This variation is used as \ref element::Element<dim>.ext_stress_coeff_. */

	// make sure that the solver doesn't have any added load nor eigenstrain
	// fields, otherwise the computed coefficients won't be the right ones
	solver.clear();

	// add the unit load and solve for the elastic fields
	solver.add_load_increment(1.0);
	solver.solve();

	// use the computed local stress as each element's ext_stress_coeff
	auto &stress = solver.get_stress();
	for(auto &element : elements)
		element->ext_stress_coeff(stress[element->number()]);

	// clear again the solver so this computation doesn't affect the future use
	// of the input solver outside this function
	solver.clear();
}


template<int dim>
void calculate_local_stress_coefficients(
	element::Vector<dim> &elements, elasticity_solver::Solver<dim> &solver)
{
	/*! Calculate, using the input elascity solver, the rank-4 tensor \f$
	 * \mathbb{S} \f$ that relates linearly local stress variations and
	 * eigenstrain increments. This tensor is set as
	 * \ref element::Element<dim>::S_. The tensor is computed individually for
	 * each element in the input vector. */

	// See the documentation of mepls::utils::compute_voigts_stiffness for a
	// detailed description of the method. Essentially, we compute the
	// components of S by solving a linear system of equations where the
	// componets of S are the unkowns, the local stress variations are the
	// left-hand side and the local eigenstrain increments form the coefficients
	// of the matrix.

	// make sure that the solver doesn't have any added load nor eigenstrain
	// fields, otherwise the computed coefficients won't be the right ones
	solver.clear();

	// create 3 independent eigenstrain increments: [[1,0],[0,0]], [[0,0],[0,1]]
	// and [[0,1],[1,0]]
	dealii::SymmetricTensor<2, dim> eigenstrain_a;
	eigenstrain_a[0][0] = 1.;

	dealii::SymmetricTensor<2, dim> eigenstrain_b;
	eigenstrain_b[1][1] = 1.;

	dealii::SymmetricTensor<2, dim> eigenstrain_c;
	eigenstrain_c[0][1] = 1.;

	// We solve the linear system of equations using deal.II tools. For that,
	// we allocate the necessary objects that will be used repeatedly

	dealii::Vector<double> f(9); // rhs (stresses)
	dealii::Vector<double> u(9); // unkowns (S compnents)
	dealii::FullMatrix<double> M(9, 9); // matrix (eigenstrains)

	// we want S for each element
	for(auto &element : elements)
	{
		unsigned int n = element->number();

		M_Assert(n < elements.size(), "");

		// solve the 3 elastic problems
		solver.add_eigenstrain(n, eigenstrain_a);
		solver.solve();
		auto stress_a = solver.get_stress()[n];
		solver.clear();

		solver.add_eigenstrain(n, eigenstrain_c);
		solver.solve();
		auto stress_c = solver.get_stress()[n];
		solver.clear();

		solver.add_eigenstrain(n, eigenstrain_b);
		solver.solve();
		auto stress_b = solver.get_stress()[n];
		solver.clear();

		// We solve the linear system of equations to obtain S. We use a
		// function that computes the rank-4 stiffness tensor. This is ok,
		// assuming that S has the same symmetries.
		auto S_voigt = utils::tensor::compute_voigts_stiffness<dim>(stress_a, stress_b, stress_c,
																	eigenstrain_a, eigenstrain_b,
																	eigenstrain_c, f, u, M);

		AssertIsFinite(S_voigt[0][0]);
		AssertIsFinite(S_voigt[0][1]);
		AssertIsFinite(S_voigt[0][2]);
		AssertIsFinite(S_voigt[1][1]);
		AssertIsFinite(S_voigt[1][2]);
		AssertIsFinite(S_voigt[2][2]);M_Assert(not std::isnan(S_voigt[0][0]), "");M_Assert(
			not std::isnan(S_voigt[0][1]), "");M_Assert(not std::isnan(S_voigt[0][2]), "");M_Assert(
			not std::isnan(S_voigt[1][1]), "");M_Assert(not std::isnan(S_voigt[1][2]), "");M_Assert(
			not std::isnan(S_voigt[2][2]), "");

		// write S as rank-4 and set it to the element
		element->S(utils::tensor::voigt_to_standard_rank4<dim>(S_voigt));
	}

}


template<int dim>
void calculate_local_stress_coefficients_central(
	element::Vector<dim> &elements, elasticity_solver::Solver<dim> &solver)
{
	/*! Calculate, using the input elascity solver, the rank-4 tensor \f$
	 * \mathbb{S} \f$ that relates linearly local stress variations and
	 * eigenstrain increments. This tensor is set as
	 * \ref element::Element<dim>::S_. The tensor is computed once for the element
	 * in the center of the domain and then reused for all the other elements.
	 * In this way, we avoid the expensive process of computing the tensor for
	 * every element. This is possible when the system has elastic homogeneous
	 * properties.
	 *
	 * @note the presence of surfaces breaks the material homogeneity since the
	 * elastic response near the surfaces differs from the bulk. If the system
	 * is small and does not have periodic boundaries, the presence of the
	 * surfaces might become important. In this case, it is recommended to use
	 * \ref calculate_local_stress_coefficients instead. */

	// NOTE see calculate_local_stress_coefficients(...) for further
	// documentation

	solver.clear();

	// compute the element number whose spatial coordinates lie in the center of
	// the system
	unsigned int Nx = solver.get_Nx();
	unsigned int Ny = solver.get_Ny();
	unsigned int cx = std::floor(Nx / 2);
	unsigned int cy = std::floor(Ny / 2);
	unsigned int n = cy * Nx + cx;

	dealii::SymmetricTensor<2, dim> eigenstrain_a;
	eigenstrain_a[0][0] = 1.;

	dealii::SymmetricTensor<2, dim> eigenstrain_b;
	eigenstrain_b[1][1] = 1.;

	dealii::SymmetricTensor<2, dim> eigenstrain_c;
	eigenstrain_c[0][1] = 1.;

	dealii::Vector<double> f(9);
	dealii::Vector<double> u(9);
	dealii::FullMatrix<double> M(9, 9);

	solver.add_eigenstrain(n, eigenstrain_a);
	solver.solve();
	auto stress_a = solver.get_stress()[n];
	solver.clear();

	solver.add_eigenstrain(n, eigenstrain_c);
	solver.solve();
	auto stress_c = solver.get_stress()[n];
	solver.clear();

	solver.add_eigenstrain(n, eigenstrain_b);
	solver.solve();
	auto stress_b = solver.get_stress()[n];
	solver.clear();


	auto S_voigt = utils::tensor::compute_voigts_stiffness<dim>(stress_a, stress_b, stress_c,
																eigenstrain_a, eigenstrain_b,
																eigenstrain_c, f, u, M);

	AssertIsFinite(S_voigt[0][0]);
	AssertIsFinite(S_voigt[0][1]);
	AssertIsFinite(S_voigt[0][2]);
	AssertIsFinite(S_voigt[1][1]);
	AssertIsFinite(S_voigt[1][2]);
	AssertIsFinite(S_voigt[2][2]);M_Assert(not std::isnan(S_voigt[0][0]), "");M_Assert(
		not std::isnan(S_voigt[0][1]), "");M_Assert(not std::isnan(S_voigt[0][2]), "");M_Assert(
		not std::isnan(S_voigt[1][1]), "");M_Assert(not std::isnan(S_voigt[1][2]), "");M_Assert(
		not std::isnan(S_voigt[2][2]), "");

	// use the computed S for all the elements
	for(auto &element : elements)
		element->S(utils::tensor::voigt_to_standard_rank4<dim>(S_voigt));

}


/*! This struct is used for output purposes only. It dumps a complex 
 * \ref element::Element object into a simple struct of scalar values, which
 * contains the configuration with which a certain element was created. In this
 * way, it can be easily written into, e.g., hdf5 datasets. */
template<int dim>
struct SetupRow;

/*! Instantiation of ElementSetupRow for dimension 2. */
template<>
struct SetupRow<2>
{
	unsigned int type;
	/*!< Element type. See \ref element::Element<dim>.type_. */

	float prestress_00;
	/*!< Component 00 of \ref element::Element<dim>.prestress_. */

	float prestress_11;
	/*!< Component 11 of \ref element::Element<dim>.prestress_ . */

	float prestress_01;
	/*!< Component 01 of \ref element::Element<dim>.prestress_. */

	float ext_stress_coeff_00;
	/*!< Component 00 of \ref element::Element<dim>.ext_stress_coeff_. */

	float ext_stress_coeff_11;
	/*!< Component 11 of \ref element::Element<dim>.ext_stress_coeff_. */

	float ext_stress_coeff_01;
	/*!< Component 01 of \ref element::Element<dim>.ext_stress_coeff_. */

	float C_00;
	/*!< Component 00 of \ref element::Element<dim>.C_. */

	float C_01;
	/*!< Component 01 of \ref element::Element<dim>.C_. */

	float C_02;
	/*!< Component 02 of \ref element::Element<dim>.C_. */

	float C_11;
	/*!< Component 11 of \ref element::Element<dim>.C_. */

	float C_12;
	/*!< Component 12 of \ref element::Element<dim>.C_. */

	float C_22;
	/*!< Component 22 of \ref element::Element<dim>.C_. */

	unsigned int number;
	/*!< Element number. */
};


} // namespace element

} // namespace mepls

#endif //RUN_SIM_ELEMENTS_H
