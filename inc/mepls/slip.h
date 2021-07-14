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

#ifndef RUN_SIM_SLIP_H
#define RUN_SIM_SLIP_H

#include <deal.II/base/symmetric_tensor.h>

#include <mepls/utils.h>
#include <mepls/history.h>

namespace mepls
{

namespace element
{
template<int dim>
class Element;
}


/*! This namespace contains the abstract slip class @ref mepls::slip::Slip<dim>, whcih provides a
 * common interface for derived slip classes. */
namespace slip
{

/*! Slip objects define the physics of activation of the slip systems within
 * mesoscale elements. When a slip system becomes active, a local slip event takes place
 * (see \cite DFCastellanos_CRP). This abstract class provides a common interface for derived
 * slip classes, which implement specific activation rules.
 * Essentially, slip systems are defined by a slip plane, an effective shear
 * stress on that plane, and the critical effective shear stress known as the
 * threshold that defines when the system becomes mechanically activated. The
 * interplay between these quantities is implemented by derived classes. This class
 * serves as a common interface for the slip classes.
 *
 * @note this class defines only the slip activation rules but **does not**
 * represent an active slip event, which is represented by @ref mepls::event::Plastic<dim>.
 */
template<int dim>
class Slip
{
  public:

	Slip()
	{
		/*! Constructor. */
	}

	Slip(Slip *original)
	{
		/*! Copy constructor. Copy the value of all the data members except the
		 * pointer to the parent element. Since the copied slip will in general
		 * be added to a different elemetn, this pointer should be set by that
		 * element when the slip is added to it. */

		eff_shear_stress = original->eff_shear_stress;
		pressure = original->pressure;
		threshold = original->threshold;
		barrier = original->barrier;
		activation_rate = original->activation_rate;
		angle = original->angle;
	}

	virtual ~Slip()
	{
		/*! Virtual destructor. */
	}

	virtual void update() = 0;
	/*!< Update the state of the slip object using on the state of the parent
	 * mesoscale element, @ref parent. */

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

		return make_copy_impl();
	}

	double eff_shear_stress = 0.;
	/*!< Effective shear stress \f$ \tau \f$ in the slip system. It is
	 * responsible for lowering the slip activation barrier, @ref barrier. */

	double pressure = 0.;
	/*!< The pressure. Since mesoscale elements have a single stress tensor, all the
	 * slip systems owned by the same element share the same pressure. */

	double threshold = 0.;
	/*!< Local yield threshold \f$ \tau^{\rm c} \f$. It is the critical resolved shear stress
	 * \f$ \tau \f$ that the slip system needs to reach to become active (mechanical unstable). */

	double barrier = 0.;
	/*!< Stress distrance to threshold, \f$ \Delta\tau^{\rm c} = \tau^{\rm c} - \tau \f$
	 * neccesary to be overcome before the slip system becomes mechanically activated. */

	double activation_rate = 0.;
	/*!< Rate \f$ \nu \f$ for thermal slip activation. */

	double angle = 0.;
	/*!< Orientation of the slip plane (in radians). */

	element::Element<dim> *parent;
	/*!< Parent mesoscale element element::Element to which the slip object
	 * belongs. When the state of the slip oject is updated via @ref update(),
	 * this element is accessed.
	 *
	 * @warning parent elements contain pointers to the owned slip objects, thus
	 * deletion of slips objects must be done carefully. Usually, allocation and
	 * deletion of slip objects is only done within the implementation of
	 * element::Element objects. */

  private:

  	virtual Slip<dim> *make_copy_impl()
	{
		/*! Dynamically allocate a new object of the derive class type and
		 * return a pointer to it. */

		M_Assert(false, "Copying of the current type of slip not implemented.");

		abort();
	}
};

} // namespace slip


} // namespace mepls

#endif //RUN_SIM_SLIP_H
