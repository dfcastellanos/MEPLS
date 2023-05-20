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

#ifndef __EVENT_H_
#define __EVENT_H_

#include <deal.II/base/symmetric_tensor.h>


namespace mepls
{

namespace slip
{
template<int dim>
class Slip;
}

/*! @namespace mepls::event
 * @brief This namespace contains event classes, which represent
 * the occurrence of different kinds of discrete events within the system. */
namespace event
{

/*! @class mepls::event::Plastic
 * @brief This class represents a local plastic deformation event. This plastic event is the
 * continuum-mechanics representation of a localized slip event in the material's microstructure. */
template<int dim>
struct Plastic
{
	Plastic(slip::Slip<dim> *slip_)
		:
		slip(slip_),
		acting_stress(slip_->parent->stress()),
		element(slip_->parent->number())
	{
		/*! Constructor.
		 * @warning the input pointer can be invalidated if the slip's parent
		 * element @ref slip::Slip::parent gets its structural properties renewed
		 * by calling to element::Element::renew_structural_properties, and the
		 * implementation of such function by derived element classes involves
		 * removing the existing slips. */
	};

	Plastic()
	{
		/*! Constructor.
		 *
		 * @warning after constructing this object, the member slip remains
		 * uninitialized. */
	};

	slip::Slip<dim> *slip;
	/*!<  Pointer to the slip object which has been activated. */

	double dplastic_strain = 0.;
	/*!< Scalar effective plastic shear strain (see @ref
	 * element::Element<dim>::integrated_vm_eigenstrain). */

	double dtime = 0.;
	/*!< Duration of the event. */

	double slip_threshold = 0.;

	dealii::SymmetricTensor<2, dim> eigenstrain;
	/*!< Eigenstrain added to the parent element (see @ref slip::Slip<dim>::parent) of
	 * the activated slip object due to the plastic event. See @ref
	 * element::Element.eigenstrain. */

	dealii::SymmetricTensor<2, dim> acting_stress;
	/*!< Stress tensor acting on the parent element (see @ref slip::Slip<dim>::parent)
	 * at the moment of the platic event. */

	unsigned int activation_protocol = 0;
	/*!< Dynamics protocol (see @ref dynamics::Protocol) by which the plastic
	 * event has been triggered. */

	unsigned int element = 0;
	/*!< Number of the element in which the plastic deformation takes place. */

	bool renew_slip_properties = true;
	/*!< Renew the slip systems owned by the active slip's parent element. */

	bool renew_elastic_properties = false;
	/*!< Renew the elastic properties of the active slip's parent element. */
};


/*! @class mepls::event::Driving
 * @brief This class represents changes in the external driving conditions. */
template<int dim>
struct Driving
{
	Driving()
		:
		dext_stress(0.),
		dpressure(0.),
		dtotal_strain(0.),
		dload(0.),
		dtime(0.),
		activation_protocol(0)
	{
		/*! Constructor. Initialize members to 0. */
	};

	double dpressure;

	double dext_stress;
	/*!< Change of the applied external stress. */

	double dtotal_strain;
	/*!< Change of the total strain. */

	double dload;
	/*!< Change of the external load. */

	double dtime;
	/*!< Time interval over which the change of driving conditions occurs. */

	unsigned int activation_protocol;
	/*!< Dynamics protocol (see @ref dynamics::Protocol) responsible of
	 * varying the external driving conditions.*/
};


/*! @class mepls::event::RenewSlip
 * @brief This class represents changes in the local microstructural properties. */
template<int dim>
struct RenewSlip
{
	RenewSlip(slip::Slip<dim> *slip)
		:
		element_number(slip->parent->number()),
		threshold(slip->threshold),
		slip_angle(slip->angle)
	{
		/*! Constructor. Initialize the pointers. */
	};

	unsigned int element_number;
	/*!< Number of the parent element of the slip object. */

	double threshold;
	/*!< Copy of @ref mepls::slip::Slip<dim>::threshold. */

	double slip_angle;
	/*!< Copy of @ref mepls::slip::Slip<dim>::angle. */
};

} // namespace event


} // namespace mepls

#endif //__EVENT_H_
