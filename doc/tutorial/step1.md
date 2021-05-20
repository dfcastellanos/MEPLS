@page Step1

# Step 1

[< previous ](@ref StepNULL) | [next >](@ref Step2)

## Introduction

In this tutorial step, we will introduce the element and the slip modules. These modules contain 
the classes @ref mepls::element::Element<dim> and @ref mepls::slip::Slip<dim>, which are 
the building blocks of any model implemented with MEPLS. The element class represents the 
mesocale subdomains
in which we discretize the simulated material. The slip class represents the internal mechanisms 
present within the elements by which slip events occur. First, you will learn the basics and 
explore the essential parts of their interfaces. Next, you can find the example program 
commented step by step. After that, you can find the results of running the 
program. (in this case, the results correspond just to the information printed on the 
screen). Finally, you can find the entire program.


The class @ref mepls::element::Element<dim> is abstract and defines an interface common to 
all the element classes. Some of its virtual functions are impelemented by derived element objects
with the goal of simulating materials with different microstructural properties. For this tutorial, 
we will create an object of class @ref example::element::Scalar<dim>. Since we rely on
the standard element interface, the example also applies to any other element classes. 


Each element posses one or more slip objects of a class derived from the @ref 
mepls::slip::Slip<dim>, which provides an abstract interface to all slip classes. Slip objects 
represent the mechanisms by which a plastic deformation takes place in a specific type of material.
In other words, they define how a slip event occurs within a mesoscale element. It's important to
 remark that they represent the "how." but not the "when," which is controlled elsewhere by the dynamic 
protocols that we will see in later tutorial steps. Since the elements represent material 
subdomains, these slip events consequently represent plastic deformation which 
is spatially localized. 

The example class @ref example::element::Scalar<dim> implements 2D single-slip crystal plasticity, 
where slip events can occur along a single plane of fixed orientation. This behavior is 
implemented by slip objects of class @ref example::slip::Scalar<dim>, which implements some of the 
virtual functions of the base class @ref mepls::slip::Slip<dim>. In this case, we can interpret 
the slip objects as slip systems. However, in general, this need not be so. For example, a slip 
object can implement J2 plasticity, where no specific crystallographic plane or predefined slip angle 
exists. Nonetheless, for simplicity, we will often talk about slip systems regardless of the 
actual physical mechanisms.

For a certain tensorial stress \f$ \boldsymbol{\Sigma} \f$, the resolved shear stress \f$ 
\tau \f$ on the slip plane is

\f[ 
    \tau = M: \boldsymbol{\Sigma}
\f]

where the tensor \f$ M \f$ is given by \f$ M = ( \boldsymbol{s} \otimes \boldsymbol{n} +  
\boldsymbol{n} \otimes \boldsymbol{s})/2\f$. Here, \f$ \boldsymbol{n} \f$ is the unit vector to the slip 
plane and \f$ \boldsymbol{s} \f$ is the unit vector of the slip direction. In 2D, we can write
this tensor as a function of the slip angle \f$ \theta \f$ (also the angle of maximum shear) as

\f[
M(\theta) = \frac{1}{2}\begin{pmatrix} -\textrm{sin}2\theta & 
\textrm{cos}2\theta \\ \textrm{cos}2\theta & \textrm{sin}2\theta \end{pmatrix}
\f]

A tensorial plastic strain increment 
\f$ \Delta\boldsymbol{\varepsilon}_{\textrm pl} \f$ associated with a slip event within that 
plane is given by

\f[ 
    \Delta\boldsymbol{\varepsilon}_{\textrm pl} = \gamma M(\theta)
\f]

where \f$ \gamma \f$ is the shear strain increment. 

The class @ref example::element::Scalar<dim> defines slips in two direction, a positive one
associated with a shear angle \f$ \theta=0 \f$ and its reversed, associated with \f$ \theta=\pi/2
 \f$
 
 
## The commented program

First, we include the headers that are necessary for this tutorial step. The 
element and slip classes that we will use are defined in the header 
example.h. The header mepls/utils.h contains utility functions
and classes that come in handy many times. The header `random` belongs to the 
C++ standard library, and we will use it to import the random number generator. 
Afterwards, we create the test element that we will use throughout this tutorial.
 

```cpp
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
```

To access and modify the state of the element, we can set and get the value of 
its members as follows: `element.X()` returns the value of member `X`, while `element.X(X_)` 
sets the value `X_` to it. For example, to set the stress tensor of the element, first we create
it, specify the value of its components and then we pass it to the element. The same for the 
elastic properties:

```cpp
   // A dealII's rank-2 symmsetric tensor. In this case, the stress tensor dealii::SymmetricTensor<2,dim> stress;
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
```

In the case of @ref example::element::Scalar<dim>, elements have a single slip plane with 
a fixed orientation of @ref example::slip::Scalar<dim>::angle radians. MEPLS associates a 
different slip object to each possible slip angle. Thus, to consider the possibility of a slip 
event in the opposite direction but within the same plane, the element also has a second slip 
object associated with the slip angle \f$ \pi/2 \f$.  
    
Let's take a look: 

```cpp
   //  We can iterate over the slips present in the element std::cout << "Slip angles = ";
   for(auto &slip : element)
      std::cout << slip->angle << ", ";
   std::cout << "\n\n";

   // We can ask for the number of slip objects
   std::cout << "Number of slip systems = " << element.size() << "\n\n";

   // and we can get a pointer to a specific slip as
   unsigned int slip_number = 1;
   mepls::slip::Slip<dim> * slip = element.slip( slip_number );
```

The effective shear stress, @ref mepls::slip::Slip<dim>::eff_shear_stress represents the driving
force responsible of activating slip events. In the case of crystal plasticity, this corresponds
to the resolved shear stress on the slip plane.

Slip systems have a critical shear stress called threshold, @ref mepls::slip::Slip<dim>::threshold.
For @ref example::slip::Scalar<dim>, the threshold is randomly distributed 
to represent microstructural disorder. In the case of crystals disorder is originated by 
imperfections, such as lattice dislocations. We remark that different classes of elements could 
include other sources of randomness, such as randomly oriented local slip planes or a varying 
number of them, heterogeneous elastic properties, local damage accumulation etc.  

Another important member of the slip class is the barrier, @ref mepls::slip::Slip<dim>::barrier. The
barrier measures the difference between the threshold and the effective shear stress. It plays a
role similar to a local yield function.


```cpp
    std::cout << "Properties of slip system #" << slip_number << ": ";
    std::cout << "\n    * angle = " << slip->angle
              << "\n    * threshold = " << slip->threshold;

   //  The slips have access to their parent's state, therefore the shear stress in the slip
   //  plane is already set
   std::cout << "\n    * shear stress = " << slip->eff_shear_stress;
   std::cout << "\n    * barrier = " << slip->barrier << "\n\n";
```


We have discussed above the role played by the structural disorder. However, a 
fundamental part of mesoscale models is the structural evolution, i.e., how the 
disorder evolves. The continuum mechanics fields present in the system are element-wise 
homogeneous. By doing so, the mesoscale description suffers a loss of information with respect to
 the underlying microscale truth. Such loss means that we have imperfect knowledge about the 
 elements' state and internal dynamics, leading to statistical uncertainty. To capture such 
uncertainty, structural properties evolve according to stochastic processes. 
The reason behind that evolution is plastic deformation, which leads to 
permanent atomistic rearrangements. Since the local atomistic configuration is 
at the origin of the local microstructural properties, it is natural to 
consider that plasticity modifies the microstructure. In the case, e.g., of 
crystal lattices, structural evolution is associated to changes in the 
configuration of a dislocation network. 

In MEPLS, we can inform an element that its structural properties must be 
renewed by calling its function @ref mepls::element::Element<dim>::renew_structural_properties():

```cpp
   element.renew_structural_properties();
```

Which properties are renewed, and in which specific way they are renewed, 
depends on the specific element class. In this case, @ref example::element::Scalar<dim> implements
 are very simple renewal process: the slip threshold is renewed by drawing a new value from its 
 probability distribution in an independent manner. Since the slip planes do not change, the slip
 angles remain fixed.  

To see this, we can recheck the slip properties. However, we must note that 
slip objects are, in general, not permanent, i.e., their parent element might 
delete them and create new ones (this is up to the specific element class to 
decide how to handle slip objects). In that case, existing pointers become 
invalidated. Therefore, we create a new pointer to one of the slips before accessing 
its members:


```cpp
   // the previous slip doesn't exist anymore
   slip = element. slip( slip_number );

   std::cout << "Properties of the *new* slip system #" << slip_number << ": "
            << "\n    * angle = " << slip->angle
            << "\n    * threshold = " << slip->threshold
            << "\n    * shear stress = " << slip->eff_shear_stress
            << "\n    * barrier = " << slip->barrier << "\n\n";
}
```

As expected, since we didn't change the parent's stress, the `eff_shear_stress` 
remains the same, however, the threshold and the barrier are different.

In the netx tutorial step (@ref Step2), we will keep exploring the interfaces of the element and 
slip classes. We will see e.g. how to get the plastic strain increments associated to a slip 
event and how to add it to an element.



## Results

This is the terminal output:

```
Stress = 1 0.5 0.5 -1

Elastic energy = 0.625

Slip angles = 0, 1.5708, 

Number of slip systems = 2

Properties of slip system #1: 
    * angle = 1.5708
    * threshold = 0.874838
    * shear stress = -0.5
    * barrier = 1.37484

Properties of the *new* slip system #1: 
    * angle = 1.5708
    * threshold = 1.00718
    * shear stress = -0.5
    * barrier = 1.50718
```



## The full program

```cpp
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
```

[deal.II]: https://www.dealii.org/