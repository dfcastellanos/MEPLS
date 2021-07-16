@page Step2

# Step 2

[< previous ](@ref Step1) | [next >](@ref Step3)

### Table of contents

- [Introducion](#introducion_2) 
- [The commented program](#comented_program_2)
    - [The solver](#the_solver)
    - [The system](#the_system)
    - [Plastic events](#plastic_events)
    - [Driving events](#driving_events)
    - [Cleaning up](#cleaning)
- [Results](#results_2)
- [The complete program](#full_2)

# Introduction{#introducion_2}

In this tutorial step, we will introduce the solver and system modules. These modules contain the 
classes @ref mepls::elasticity_solver::Solver<dim> and @ref mepls::system::System<dim>, 
whose goal is to make the mesoscale elements work together. Among other things, we will see how 
to get the plastic strain increments associated with a slip event and how to add it to an element.

The solver class calculates elastic fields using the Finite Element Method. For that, it uses a 
mesh where each finite element is associated with one mesoscale element. The system class brings 
together the solver and the elements and provides higher-level functions to manage them based on 
events. It will also automatically take care of recording the events and the evolution of many 
system properties.

We will continue using the same slip and elements classes introduced in @ref Step1, i.e.,@ref 
mepls::slip::Slip<dim> and  @ref mepls::element::Element<dim>.


# The commented program{#comented_program_2}

First, we include the headers that are necessary for this tutorial step. We will include the same
headers of the previous step (see @ref Step1) plus the headers containing the solver and system 
classes. The event classes are defined in the event header.
 

```cpp
#include <example.h>
#include <mepls/utils.h>
#include <random>

// new headers
#include <mepls/solver.h>
#include <mepls/system.h>
#include <mepls/event.h>
```

Some things will be necessary for creating the elements, as in the previous tutorial step: 

```cpp
int main()
{
   constexpr unsigned dim = 2;
    std::mt19937 generator(1234567);

   double G = 1.;
   double nu = 0.3;
   dealii::SymmetricTensor<4,dim> C = mepls::utils::tensor::make_isotropic_stiffness<dim>(G, nu);
```

In the previous step, @ref Step1, we showed an element's elastic behavior in isolation. We 
focus now on how it behaves when embedded in an elastic matrix formed by other elements. 

MEPLS handles the elements as a 1D vector, independently of the dimensions of the problem (
in this case, 2D). The reason is that elements behave only according to their own (local) state, and
their actual spatial location only matters when computing the elastic fields. That 
means that only the solver (internally) cares about spatial coordinates, but outside the solver
most of the MEPLS classes don't. Let's discretize the domain into a lattice of 16x16 mesoscale 
elements:

```cpp
   unsigned int Nx = 16;
   unsigned int Ny = 16;

    // this type is an alias for std::vector<mepls::element::Element<dim>>
   mepls::element::Vector<dim> elements;
   
   for(double n = 0; n < Nx * Ny; ++n)
   {
      example::element::Scalar<dim>::Config conf;
      conf.number = n;

      auto element = new example::element::Scalar<dim>(conf, generator);
      element->C( C );

      elements.push_back(element);
   }
```

## The solver{#the_solver}

The solver class calculates elastic fields using the Finite Element Method. For that, it
uses a mesh where each finite element is associated with one mesoscale element. The
average stress tensor computed over each of the finite elements is used as the stress tensor
of the associated mesoscale element. The class @ref mepls::elasticity_solver::Solver<dim> provides 
an abstract interface common to all the solver classes. Some of its virtual functions are 
implemented by derived solver classes, allowing to compute the elastic fields under different 
conditions but always relying on the standard interface provided by the base class. 

In this tutorial, we will use the solver @ref mepls::elasticity_solver::LeesEdwards<dim>,
but since it shares an interface common to all the solvers in the namespace @ref 
mepls::elasticity_solver, this example is of wide applicability.

The solver @ref mepls::elasticity_solver::LeesEdwards<dim> computes elastic fields with 
bi-periodic boundary conditions, under the action of a load that induces an external shear stress
field \f$ \boldsymbol{\Sigma}_{\rm ext} = \Sigma_{\rm xy} (\boldsymbol{e}_{\rm x} 
\otimes \boldsymbol{e}_{\rm y} + \boldsymbol{e}_{\rm y} \otimes \boldsymbol{e}_{\textrm
x}) \f$. The finite element mesh is quadrilateral and structured, and threfore the material 
domain discretization can be understood as a lattice of mesoscale elements. 

@note currently, all the MEPLS solvers work only with structured quadrilateral meshes, and solve 
the stress equilibrium equation in static conditions (i.e., no elastic waves). Moreover, the 
deformation is always computed with respect to the original pristine mesh, so the solver does
not re-mesh.

We will use a mesh of 16x16 elements to match the size of the vector of mesoscale elements
created above. Each finite element has faces with a length of 1.0, which defines the units of length in the simulation.


```cpp 
   // create the solver
   mepls::elasticity_solver::LeesEdwards<dim> solver(Nx, Ny);
   
   // set the elastic properties using the ones specified for the elements
   for(auto &element : elements)
      solver.set_elastic_properties(element->number(), element->C());
      
    // this operation builds the FEM system. After this is called, the solver is
    // ready to be used and its mesh geometry or loading model cannot be altered
   solver.setup_and_assembly();

```

Before continuing, let's take a very brief detour to show the basic members of the 
solver interface. Usually, this interface is used only by other MEPLS classes such that
the user needs not care about it. However, a basic knowledge will help understand how MEPLS works.


The solver computes elastic fields that arise due to an external load mechanism 
and an internal plastic strain field. The external load is composed of 3 things: the loading 
modes (pure shear, compression, bending, etc.), the driving mode (strain- or stress-controlled) 
and the load value. The load value is a float number that controls the intensity of the applied 
load. Its exact meaning depends on the specific solver implementation, but in general it is 
proportional to the applied strain or stress. In this example, we are using 
@ref mepls::elasticity_solver::LeesEdwards<dim>, which implements a pure shear loading with 
(default) strain-controlled conditions. In this case, the meaning of the load value here 
denotes the xy-component of the strain imposed by the loading mechanism. 

@note the specific meaning of the load value depends on the class of solver and loading 
conditions. For further details, see the documentation of specific solvers in @ref 
mepls::elasticity_solver.

After calling @ref mepls::elasticity_solver::Solver<dim>::setup_and_assembly(), the type of 
loading mode and driving mode is fixed, but the value of the load can be varyied by calling @ref
 mepls::elasticity_solver::Solver<dim>::add_load_increment().

In MEPLS, plastic strain is represented as an eigenstrain field. Eigenstrain denotes the deformation
that a material inclusion would undergo if it were in a stress-free situation, e.g., as cut out of 
the material matrix @cite MICHEL1999109 @cite Lifeng2014. As with the other continuum mechanic 
fields, the eigenstrain field is discretized and is element-wise homogeneous over the elements. 
When we want to compute the solution to the stress equilibrium equation under a certain 
eigenstrain field and a certain external load value, we call @ref mepls::elasticity_solver::Solver<dim>::solve().
This can be done only after @ref mepls::elasticity_solver::Solver<dim>::setup_and_assembly() has 
been called. After solving, we can get the stress field as @ref mepls::elasticity_solver::Solver<dim>::get_stress(). Let's see an 
example of this:


```cpp
   std::cout << "Testing the solver..." << std::endl;

   double delta_external_epsilon_xy = 0.1;
    solver.add_load_increment( delta_external_epsilon_xy );

    dealii::SymmetricTensor<2,dim> local_eigenstrain_increment;
   local_eigenstrain_increment[0][1] = 0.2;
   unsigned int element_number = 120;
   solver.add_eigenstrain(element_number, local_eigenstrain_increment);

   solver.solve();

   const std::vector<dealii::SymmetricTensor<2,dim>> & stress_field = solver.get_stress();
   const std::vector<dealii::Tensor<2,dim>> & def_grad = solver.get_deformation_gradient();

    double external_stress = solver.get_external_stress();
    double total_strain = solver.get_total_strain();

   // just tau = G * 2eps_xy
   std::cout << "    * external stress = " << external_stress << std::endl;

   // just eps_xy
   std::cout << "    * total strain = " << total_strain << std::endl;

   std::cout << "Cleaning the solver..." << std::endl;

   // this call resets the solver state, so the changes we made to load and the
   // eigenstrain field disappear
   solver.clear();    
```
At the end, we call @ref mepls::elasticity_solver::Solver<dim>::clear(), so all these 
example changes that we have just done are removed. The solver is back to its clean state, 
exactly as it was just after calling @ref mepls::elasticity_solver::Solver<dim>::setup_and_assembly().

## The system{#the_system}

After this quick example of how the solvers operate, let's retake our main goal now: building the
system of mesoscale elements.

In the previous tutorial step, we saw how to set the elements' structural properties. Now, we will
inform the elements about how the external loading mechanism changes their stress state upon 
variation of the external load value. To do it, we call the function 
@ref mepls::element::calculate_ext_stress_coefficients(). This function computes the 
local stress change induced by a load increment of value 1.0 and uses it to set the value of
@ref mepls::element::Element<dim>::ext_stress_coeff. Knowing @ref 
mepls::element::Element<dim>::ext_stress_coeff will allow us to ask to each 
specific slip systems, which is the minimum externalload increment necessary to unstabilize it, i.e.
to make its barrier 0. This is also known as the slip's critical external load increment.
The value of the critical load increment is obtained by calling 
@ref mepls::slip::Slip<dim>::get_critical_load_increment().


```cpp
   mepls::element::calculate_ext_stress_coefficients(elements, solver);

   // let's take some element
   auto element = elements[40];
   
   // let's take the first slip of that element
   auto slip = element->slip(0);

   // this is the external load increment that is necessary to make this slip's sytem
   // shear stress match its threshold, that is to make the barrier equal 0   
   double critical_incr = slip->get_critical_load_increment();
   
   std::cout << "\nSlip's system critical load increment = " << critical_incr << std::endl;
```

We can ask that specific slip system for the eigenstrain increment associated with a slip event 
on its plane. Such increment represents, at the mesoscale, the effects of a slip event 
occurring within the element. We can obtain it by calling @ref 
example::slip::Scalar<dim>::get_eigenstrain_increment().
How this eigenstrain is calculated depends on the specific class of
slip object. Many times, that eigenstrain is related to the local acting stress by some 
rule. For example, a simple such rule might define that slips amplitudes must lead to local 
shear stress drops of a 10%. Applying such a rule requires prior knowledge of how local plastic 
deformation leads to changes in local stress. In other words: if an element suffers a change in its 
eigenstrain, what local stress changes does it cause? This information can be obtained by calling @ref 
mepls::element::calculate_local_stress_coefficients(). This function applies different local 
eigenstrain changes and computes the resulting local stress changes. With these changes, the 
function calculates a rank-4 tensor relating them and sets it in the elements with @ref 
mepls::element::Element<dim>::S.

In general, this tensor is different for each element if there are elastic inhomogeneities, 
different element shapes or surfaces. In this example, the system is elastically homogeneous and
has periodic boundaries. On the other hand, the elements have the same shape and size 
(since we use a structured finite element mesh). In this case, the tensor is the same for 
all the elements, and it would be a waste of effort to compute it for each element, which can 
take considerable time. To avoid computing it everywhere, we call instead @ref 
mepls::element::calculate_local_stress_coefficients_central, which computes it in the element 
located in the center of the domain and reuses it for every other element.

@note the class of slips used in this tutorial, @ref example::slip::Scalar<dim>, implement fixed 
slip amplitudes with a value given by @ref example::slip::Scalar<dim>::gamma, so knowing @ref 
mepls::element::Element<dim>::S is not strictly neccesary. However, it is illustrative to show it
here.


```cpp
   mepls::element::calculate_local_stress_coefficients_central(elements, solver);

    // this increment is to be added to the slips parent element to represents the effects
    // of the slip event
   dealii::SymmetricTensor<2,dim> eigenstrain_incr = slip->get_eigenstrain_increment();
   
   std::cout << "Slip's system local eigenstrain increment = " << eigenstrain_incr << std::endl;
```

As expected, the eigenstrain increment has only a xy-component due to the orientation of the 
slip plane for @ref example::slip::Scalar<dim>, as explained in the previous tutorial (@ref Step1).
This eigenstrain could be added to the solver, as we saw above. 
However, if we want the elements to be aware of their plastic deformation, we should also add
the increment to the slip's  parent element with @ref 
mepls::element::Element<dim>::add_eigenstrain().
This call does not perform changes in the system, but it allows the elements to have an updated 
knowledge of their plastic deformation. This matters if we want to request, e.g., the value of 
their local von Mises plastic strain. 

@note The reason why we need to add the increment to the solver 
and the elements separately is a design decision, which allows different MEPLS modules to remain
as independent and uncoupled from each other as possible. This design also gives the users a 
great deal of freedom when implementing models since they can be very specific about how the 
different modules should interact with each other.

However, the fact that we must perform several function calls from different objects to 
represent a single operation (adding a slip event's plastic deformation) means that there's still
a lot of room for abstraction. Therefore,  we take another step in the abstraction of the problem 
and introduce the @ref mepls::system::System<dim> class, which brings together the solver and the 
elements and provides higher-level functions to manage them based on events. The events 
classes are defined in the namespace @ref mepls::event, and represent
 the occurrence of local slip events (@ref mepls::event::Plastic<dim>) and variations of the 
 external loading conditions (@ref mepls::event::Driving<dim>). The 
 system class also takes care of recording the events and the evolution of the magnitudes of 
 interest, which will be used for output purposes in the following tutorial. Again, the class @ref 
mepls::system::System<dim> provides an abstract interface common to all the system classes, while 
its derived classes implement specific ways of handling events. 

To understand what this means, let's remember that elements can renew their structural properties 
when they are asked to do so (see @ref Step1). Among other things, the system class will be 
responsible for telling the elements when to do it. For example, a certain system class might 
renew an element's properties whenever a slip event occurs within that element. However, 
different models might choose to handle this differently. For example, another model might 
require that slip events induce a renew of the element's neighborhood too. Or a model might 
require that a certain plastic strain increment must not be added not to one element, but 
distributed over the neighbohood. The possibilities are vast, and the system class is the place 
where this behavior is implemented. In summary, the system class offers many possibilities for 
defining what happens to the elements when event occur, and how and when the solver should update
 the elastic fields.

In this tutorial, we will use the system class @ref mepls::system::Standard<dim>. This class 
provides a simple behavior: when a slip event occurs within an element, that element's 
structural properties are renewed, and a plastic strain increment is added to that element. 
Then, the elastic fields are re-computed, taking into account the changes in the plastic strain 
field and loading conditions (the loading conditions change after plastic deformation because, 
under strain-controlled conditions, plastic deformation leads to drops of the external stress).
 
Let's create the system, and then get a handle to an internal object of class @ref 
mepls::system::MacroState<dim> that keeps track of the 
macroscale state of the system (that is, of the global properties such as external stress, total 
strain, elapsed time etc.):

```cpp  
   mepls::system::Standard<dim> system(elements, solver, generator);
   
   auto & macrostate = system.macrostate;
```

@note the system object stores references to the vector of elements, the solver and the 
generator, therefore these objects are expected to live during the entire life of the system object.

## Plastic events{#plastic_events}

Now we can easily perform slip events, and the system will take care of everything for us. To do 
it, first we create a event of class @ref mepls::event::Plastic<dim>. The event takes 
as argument the slip that we want to activate. It is performed by calling @ref 
mepls::system::System<dim>::add().

```cpp  
   mepls::event::Plastic<dim> plastic_event( slip );
   system.add(plastic_event);
```

To understand in detail what happened, let's summarize the main operations 
performed by the add function, each of which we already saw before:

1. Compute the local eigenstrain increment by calling @ref mepls::slip::Slip<dim>::get_eigenstrain_increment()
   from the active slip

2. Add that increment to the solver

3. Add that increment to the slip's parent element

4. Solve the FEM problem

5. Update the stress tensor of all the elements

6. Renew the structural properties of the parent element (note: this operation deletes the active
    slip system, invalidating the pointers to it from now on)

7. Update macrostate: external stress, total strain, plastic strain etc.

8. Record the slip event and the macrostate in the system's internal history. Record a driving 
    event (see next) associated with external stress drops or increments due to plastic deformation

These operations are implemented in the derived system class, i.e., in @ref 
mepls::system::Standard<dim>::add(mepls::event::Plastic<dim>&). After we have added the event to 
the system, let's check the effects of these operations by chekcing the parent's element and the macroscate: 

```cpp  
    std::cout << "\nAfter adding the slip event: "
              << "\n    * Local eigenstrain = " << element->eigenstrain()
              << "\n    * Local von Mises eigenstrain = " << element->integrated_vm_eigenstrain()
              << "\n    * Local stress = " << element->stress()
              << "\n    * External stress = " << macrostate["ext_stress"]
              << "\n    * Total strain = " << macrostate["total_strain"]
              << "\n    * Global von Mises plastic strain = " << macrostate["av_vm_plastic_strain"]
              << std::endl;

   system.solver.write_vtu("slip_event.vtu");              
```

We can see that the parent's element state has been updated with the same eigenstrain increment 
that we calculated above when calling @ref example::slip::Scalar<dim>::get_eigenstrain_increment() 
by ourselves. Also, the elastic effects of adding such increment have been computed, and the 
elements stress state updated accordingly. The slip event also influences the macroscale 
properties: there is a non-zero global von Mises plastic deformation and, since the total strain is 
kept fixed to 0, the external stress is now slightly negative due to
the plastic deformation.

Sometimes, it is useful to see the changes introduced by the event itself. When we create 
the event object, it is only given the active slip, and the rest of its members remain with their 
default initialization. However, when the event is added to the system, its members are 
updated with more information about the event. For example, we can check the eigenstrain increment 
associated with the event with @ref mepls::event::Plastic<dim>::eigenstrain. In this case, this 
increment is already known to us, and we won't show it. We can also check the number of the 
mesoscale element to which the increment is added by checking @ref 
mepls::event::Plastic<dim>::element. The active slip is @ref mepls::event::Plastic<dim>::slip, 
however as explained in the previous tutorial step, accessing the slip object like this is 
dangerous since it might be deleted when renewing the structural properties of its parent element
 (and indeed, it is deleted in this example).

## Driving events{#driving_events}

Let's now consider a driving event, of class @ref mepls::event::Driving<dim>. This class
allows us to control the loading mechanism by performing load increments. As with 
the slip events, the system will take care of everything for us when adding a driving 
event.

```cpp  
    mepls::event::Driving<dim> driving_event;
    
    // a load increment of 0.01
    driving_event.dload = 0.01;
    
    system.add(driving_event);
```

Let's summarize the main operations performed by the add function:

1. Add the load increment to the solver and solve the FEM problem

2. Update the stress tensor of all the elements

3. Update the external stress

4. Record the driving event and the macrostate in the system's internal history

These operations are implemented in the derived system class, i.e., in @ref 
mepls::system::Standard<dim>::add(mepls::event::Driving<dim>&). Let's take a look at the local 
stress state and the external stress after having added the driving event:

```cpp  
    std::cout << "\nAfter adding the driving event: "
              << "\n    * Local eigenstrain = " << element->eigenstrain()
              << "\n    * Local von Mises eigenstrain = " << element->integrated_vm_eigenstrain()
              << "\n    * Local stress = " << element->stress()
              << "\n    * External stress = " << macrostate["ext_stress"]
              << "\n    * Total strain = " << macrostate["total_strain"]
              << "\n    * Global von Mises plastic strain = " << macrostate["av_vm_plastic_strain"]
              << std::endl;
              
   system.solver.write_vtu("slip_and_driving_events.vtu");
   
}              
```

We see that the plastic quantities did not change by the driving event. Remember that the solver
we are using here, of class @ref mepls::elasticity_solver::LeesEdwards<dim>, is working under 
strain-controlled condition, and the load refers to the globally applied shear strain \f$ 
\varepsilon_{\rm xy} \f$. Thus, the total strain has now the value set by the load increment 
0.01 that we just added and the external stress has risen according to \f$ \tau = G 2 
\varepsilon_{\rm xy}\f$.

After the driving event object has been passed to the add function, its members have been updated. 
Thus, we could look at, e.g., the external stress change induced by the event by checking 
@ref mepls::event::Driving<dim>::dext_stress. In this case, the change in the total strain, @ref 
mepls::event::Driving<dim>::dtotal_strain, equals the change in the load, @ref 
mepls::event::Driving<dim>::dload, since we are working in strain-controlled conditions.  

## Cleaning up{#cleaning}

Lastly, we delete the elements that we created. Since they were dynamically allocated, it is our 
responsibility to delete them. In this case, the program has reached its end, so we don't really 
need to do it. But you should remember this operation if, e.g., you create vector of elements or 
systems iteratively, otherwise this will cause a memory leak.

```cpp  
	// remember that we allocated dynamically using the new operator
	// when we build a model, it will be our responsibility to delete them
	for(auto &element : elements)
		delete element;
```


# Results{#results_2}

You can compile this program as explained in [How to build](@ref HowToBuild). Then, you can run it
 as follows,

```sh
$ ./run_sim
```

This is the output:

```cpp
Testing the solver...
    * external stress = 0.198437
    * total strain = 0.1
Cleaning the solver...

Slip's system critical load increment = 0.235111
Slip's system local eigenstrain increment = -0 0.025 0.025 0

After adding the slip event: 
    * Local eigenstrain = 0 0.025 0.025 0
    * Local von Mises eigenstrain = 0.05
    * Local stress = 3.71726e-19 -0.0215013 -0.0215013 1.23909e-18
    * External stress = -0.000195313
    * Total strain = 0
    * Global von Mises plastic strain = 0.000195313
slip_event.vtu

After adding the driving event: 
    * Local eigenstrain = 0 0.025 0.025 0
    * Local von Mises eigenstrain = 0.05
    * Local stress = -5.94762e-18 -0.00150134 -0.00150134 -1.98254e-17
    * External stress = 0.0198047
    * Total strain = 0.01
    * Global von Mises plastic strain = 0.000195313
slip_and_driving_events.vtu
```

These are the generated vtu files. The one on the left shows the material after the 
slip event has been performed. We can see a local shear deformation, defined by the local 
eigenstrain increment with only xy compoenent that we added. The one on the right 
shows the material after the slip and the driving events. Note that the system undergoes a global
shearing, due to the xy shear induced by the loading conditions. The color scale corresponds to 
the xy shear stress component. The displacement of the mesh has been increased by a factor of 10 to
assist the visualization.

<center><img src="step2_events.jpg" width="50%"></center>

In this tutorial, we saw how the system class manages the elements and the solver for us. For this, 
it needs to be told what to do using slip and driving events. However, when we simulate a model, we 
don't want to manually tell the system what to do. On the contrary, we want to automatize the 
creation of events using some dynamic rules that represent a specific physical process. In the following 
tutorial, we will see how to do this and also how to access the system's evolution history for 
outputting the results.


# The complete program{#full_2}

```cpp
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

// new headers
#include <mepls/solver.h>
#include <mepls/system.h>
#include <mepls/event.h>

int main()
{

   // ------ SETTING UP THE ELEMENTS ------

   constexpr unsigned dim = 2;
    std::mt19937 generator(1234567);

   double G = 1.;
   double nu = 0.3;
   dealii::SymmetricTensor<4,dim> C = mepls::utils::tensor::make_isotropic_stiffness<dim>(G, nu);

   unsigned int Nx = 16;
   unsigned int Ny = 16;

    // this type is an alias for std::vector<mepls::element::Element<dim>>
   mepls::element::Vector<dim> elements;

   for(double n = 0; n < Nx * Ny; ++n)
   {
      example::element::Scalar<dim>::Config conf;
      conf.number = n;

      auto element = new example::element::Scalar<dim>(conf, generator);
      element->C( C );

      elements.push_back(element);
   }

   // create the solver
   mepls::elasticity_solver::LeesEdwards<dim> solver(Nx, Ny);

   // set the elastic properties using the ones specified for the elements
   for(auto &element : elements)
      solver.set_elastic_properties(element->number(), element->C());

    // this operation builds the FEM system. After this is called, the solver is
    // ready to be used and its mesh geometry or loading model cannot be altered
   solver.setup_and_assembly();



   // ------ EXPLORING THE SOLVER'S INTERFACE ------

   std::cout << "Testing the solver..." << std::endl;

   double delta_external_epsilon_xy = 0.1;
    solver.add_load_increment( delta_external_epsilon_xy );

    dealii::SymmetricTensor<2,dim> local_eigenstrain_increment;
   local_eigenstrain_increment[0][1] = 0.2;
   unsigned int element_number = 120;
   solver.add_eigenstrain(element_number, local_eigenstrain_increment);

   solver.solve();

   const std::vector<dealii::SymmetricTensor<2,dim>> & stress_field = solver.get_stress();
   const std::vector<dealii::Tensor<2,dim>> & def_grad = solver.get_deformation_gradient();

    double external_stress = solver.get_external_stress();
    double total_strain = solver.get_total_strain();

   // just tau = G * 2eps_xy
   std::cout << "    * external stress = " << external_stress << std::endl;

   // just eps_xy
   std::cout << "    * total strain = " << total_strain << std::endl;

   std::cout << "Cleaning the solver..." << std::endl;

   // this call resets the solver state so the changes we made to load and the
   // eigenstrain field dissapear
   solver.clear();


   // ------ FINISHING THE ELEMENTS SETUP ------

   mepls::element::calculate_ext_stress_coefficients(elements, solver);

   // let's take some element
   auto element = elements[40];

   // let's take the first slip of that element
   auto slip = element->slip(0);

   // this is the external load increment that is necessary to make this slip's sytem
   // shear stress match its threshold, that is to make the barrier equal 0
   double critical_incr = slip->get_critical_load_increment();

   std::cout << "\nSlip's system critical load increment = " << critical_incr << std::endl;


   mepls::element::calculate_local_stress_coefficients_central(elements, solver);

    // this increment is to be added to the slips parent element to represents the effects
    // of the slip event
   dealii::SymmetricTensor<2,dim> eigenstrain_incr = slip->get_eigenstrain_increment();

   std::cout << "Slip's system local eigenstrain increment = " << eigenstrain_incr << std::endl;


   // ------ CREATING A SYSTEM ------

   mepls::system::Standard<dim> system(elements, solver, generator);

   auto & macrostate = system.macrostate;


   // ------ ADDING A SLIP EVENT ------

   mepls::event::Plastic<dim> plastic_event( slip );
   system.add(plastic_event);

   std::cout << "\nAfter adding the slip event: "
           << "\n    * Local eigenstrain = " << element->eigenstrain()
           << "\n    * Local von Mises eigenstrain = " << element->integrated_vm_eigenstrain()
           << "\n    * Local stress = " << element->stress()
           << "\n    * External stress = " << macrostate["ext_stress"]
           << "\n    * Total strain = " << macrostate["total_strain"]
           << "\n    * Global von Mises plastic strain = " << macrostate["av_vm_plastic_strain"]
           << std::endl;

   system.solver.write_vtu("slip_event.vtu");


   // ------ ADDING A DRIVING EVENT ------

    mepls::event::Driving<dim> driving_event;

    // a load increment of 0.01
    driving_event.dload = 0.01;

    system.add(driving_event);


   std::cout << "\nAfter adding the driving event: "
           << "\n    * Local eigenstrain = " << element->eigenstrain()
           << "\n    * Local von Mises eigenstrain = " << element->integrated_vm_eigenstrain()
           << "\n    * Local stress = " << element->stress()
           << "\n    * External stress = " << macrostate["ext_stress"]
           << "\n    * Total strain = " << macrostate["total_strain"]
           << "\n    * Global von Mises plastic strain = " << macrostate["av_vm_plastic_strain"]
           << std::endl;

   system.solver.write_vtu("slip_and_driving_events.vtu");
   
    // remember that we allocated dynamically using the new operator
    // when we build a model, it will be our responsibility to delete them
    for(auto &element : elements)
        delete element;   
}
```

[deal.II]: https://www.dealii.org/
