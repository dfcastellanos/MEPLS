@page Step3

# Step 3

[< previous ](@ref Step2) | [next >](@ref Step4)

### Table of contents

- [Introducion](#introducion) 
    - [The model](#the_model)
- [The commented program](#comented_program)
    - [Setting up the system](#setting_up)
    - [The history class](#history)
    - [The stopping condition](#stopping)
    - [Rising the external load](#rising_load)
    - [Relaxing the unstable slip systems](#relaxin_slips)
    - [Writing output data](#writing_data)
- [Results](#results)
- [The complete program](#full)


# Introduction{#introducion}

This tutorial step will use the tools introduced in the previous tutorial steps to create a 
simple model. We will control the system by automatically creating events according to dynamical
rules representing the driving conditions. In the end, we will save the simulation data to an output
file and plot the results using a Python script.


### The model{#the_model}

The model we consider represents a material being driven at a very low strain rate and
temperature, the so-called athermal quasistatic limit. In this limit, the stress fields can rise
 continuously until a plastic event is triggered, but during the event, the stress does not 
rise further because the event is considered a much faster process. To understand this, we must 
take into account the characteristic duration of a plastic event, \f$ \Delta t_{\rm pl} \f$. If 
we apply an external shear strain rate \f$ \dot{\gamma}_{\rm ext} \f$, the shear stress rise 
during a plastic event is \f$ \Delta \tau = G \dot{\gamma}_{\rm ext} \Delta t_{\rm pl} \f$. 
This stress variation \f$ \Delta \tau \f$ must be compared with a typical slip threshold value given, e.g.,
 by the parameter \f$ \lambda \f$ (remember that we are using Weibull-distributed slip thresholds with scale 
\f$ \lambda \f$ and exponent \f$ k \f$, see @ref Step1). We want that, during a plastic event, the stress
does not rise significantly compared to the threshold. In this case, we obtain the following condition
for the quasistatic limit,

\f[
 \dot{\gamma}_{\rm ext} << \frac{\lambda}{G \Delta t_{\rm pl}}
\f] 


On the other hand, working in athermal conditions means that stress fluctuations induced by 
temperature are negligible compared to the slip thresholds, i.e., temperature cannot trigger 
plastic events. Energy fluctuations relate to temperature as \f$ \Delta E = k_{\rm B} T \f$ and to
stress as \f$ \Delta E \approx V_{\rm a} \Delta \tau \f$. The quantity \f$ V_{\rm a} \f$ is the 
so-called activation volume, which is the product of the typical local strain induced by a 
plastic event and the volume of the region occupied by the event. Therefore, the athermal 
limit holds if 

\f[
 T << \frac{V_{\rm a}\lambda}{k_{\rm B}}
\f] 

In our model, this limit is guaranteed because we won't include any temperature-dependent mechanism.

A significant simplification in the implementation of the model is possible by combining the 
quasistatic and the athermal limits. Since in the athermal limit plastic events cannot be 
thermally activated, the only mechanism for their activation is that local stresses overcome 
local slip thresholds. More specifically, in the case when the system is deforming elastically by
the action of the driving mechanism, at some point, a plastic (slip) event will be somewhere 
triggered by it. To simulate the activation of events in the quasistatic limit, it is enough to 
make sure that (1) the external strain increments trigger a single slip event at a time and (2) that
once a event is triggered we wait until the event (and possibly other events induced by the 
first one) ends before rising the strain again.

Therefore, to drive a system in the athermal quasistatic limit we do not need to take into 
account the time scales explictely, but just to respect their hierarchy. For this reason, we 
won't care what are the specific values of \f$ \dot{\gamma}_{\rm ext} \f$ and \f$ \Delta t_{\rm 
pl} \f$. Instead, we will care that the discrete load increments that we apply (as seen in @ref 
Step2) do not trigger many events at once (ideally, they trigger exactly 1). What is the 
criterion for the quasistatic limit in this formulation without explicit time scales? In this case,
 we require that the discrete external strain increments \f$ \Delta \gamma_{\rm ext} \f$ induce a 
stress \f$ \Delta \tau = G \Delta \gamma_{\rm ext} \f$ much smaller than \f$ \lambda \f$, i.e.,

\f[
 \Delta \gamma_{\rm ext} << \frac{\lambda}{G}
\f] 


We can summarize the driving protocol as follows: we perform discrete external strain increments 
of \f$ \Delta \gamma_{\rm ext} \f$. When a slip system becomes unstable, we keep \f$ \gamma_{\rm 
ext} \f$ fixed and perform the event (and possible, further events induced by this one). When no 
slip system is unstable, we increase the external strain again. We repeat until reaching an 
external strain limit.



# The commented program{#comented_program}
   
First, we include the headers that are necessary for this tutorial. We will include the same
headers of the previous tutorial (see @ref Step2) plus the standard library header `fstream` to 
write
the final data to an output file.
   
```cpp
#include <example.h>
#include <mepls/utils.h>
#include <mepls/solver.h>
#include <mepls/system.h>
#include <mepls/history.h>
#include <random>

// new header
#include <fstream>
```

### Setting up the system{#setting_up}

Now, we set up the elements, the solver, and the system as we saw in the previous tutorial (see @ref 
Step2). However, now we will mind the values of the simulation parameters since they will affect 
the final simulation output.

```cpp
int main()
{
   constexpr unsigned dim = 2;
    std::mt19937 generator(1234567);

   // let's consider a material with a shear modulus of 30 GPa,
   double G = 30.;
   double nu = 0.3;
   dealii::SymmetricTensor<4,dim> C = mepls::utils::tensor::make_isotropic_stiffness<dim>(G, nu);

   unsigned int Nx = 16;
   unsigned int Ny = 16;

   mepls::element::Vector<dim> elements;
   for(double n = 0; n < Nx * Ny; ++n)
   {
      example::element::Scalar<dim>::Config conf;
      conf.number = n;

      // plastic increments induce a local shear deformation of 5%
      conf.gamma = 0.05;

      // the slip thresholds have a scale of 1 GPa (see Weibull dist. for the relation between
      // lambda, k and the average value)
      conf.lambda = 1;

      // the slip thresholds have some disorder. For that, we use a not too-high k (for k->inf
      // the disorder vanishes, while it increases for k->0)
      conf.k = 6;

      auto element = new example::element::Scalar<dim>(conf, generator);
      element->C( C );

      elements.push_back(element);
   }

   mepls::elasticity_solver::LeesEdwards<dim> solver(Nx, Ny);
   for(auto &element : elements)
      solver.set_elastic_properties(element->number(), element->C());
   solver.setup_and_assembly();

   mepls::element::calculate_ext_stress_coefficients(elements, solver);
   mepls::element::calculate_local_stress_coefficients_central(elements, solver);

   mepls::system::Standard<dim> system(elements, solver, generator);
```

### The history class {#history}

To have access to the history of driving and slip events, as well as to the evolution of
different macroscale properties, we create an object of class @ref 
mepls::history::History<dim>. The macroscale properties refer to global, system-scale quantities 
such as the applied strain, the external stress, the global plastic deformation, etc. The history 
object must be passed to the system to 
inform the history about the added events. The history object will store the data in an output-friendly way, and
 we will use it at the end to write some data to a CSV file.

```cpp
   // we create a history object to record the evolution of the sytem
   mepls::history::History<dim> sim_history("Simulation_history");

   // the system needs to know the history object so that it can record the driving
   // and slip events that we add
   system.set_history(sim_history);
```

### The stopping condition {#stopping}

Now, we start the main simulation loop, in which the evolution of the system takes place. We will 
simulate until the externally applied strain reaches a target value of 5%. In each iteration, we 
will print the value of the applied strain and the external stress, whose values can be accessed 
using the `macrostate` member of the system:

```cpp
   // main simulation loop till 5% applied shear strain (epsilon_xy, not gamma)
   while( system.macrostate["total_strain"] < 0.05 )
   {

      // print some output
      std::cout << system.macrostate["total_strain"] << " " <<  system.macrostate["ext_stress"] << std::endl;
```

### Rising the external load {#rising_load}

In every main-loop iteration, we perform an increment of the applied load. 
As seen in @ref Step2, for the solver @ref mepls::elasticity_solver::LeesEdwards<dim> that we are using,
the load value is the xy-component of the the applied strain, and the system is driven under strain-controlled 
conditions. Threfore, following from the discussion in the introduction, the load increment must fulfill
\f$  \Delta \gamma_{\rm ext} << \lambda / G = 3 \cdot 10^{-2} \f$. To be on the safe side, we will 
use \f$  \Delta \gamma_{\rm ext} = 10^{-4} \f$.

```cpp
        // the applied shear strain is (instantaneously) increased by a small
        // ammout of 0.01%, and is kept fixed afterward
        mepls::event::Driving<dim> driving_event;
        driving_event.dload = 0.0001;
        system.add(driving_event);
```

### Relaxing the unstable slip systems {#relaxin_slips}
      
After adding the load increment event, the shear stress field has increased everywhere, and some 
slip systems might have become active. Since the increment was very small and respected the 
criterion for the quasistatic limit, it is very unlikely that more than one slip system is 
unstable, although the probability is not zero. As discussed in the introduction, after 
triggering an event from applying a small load increment, we must let the system locally relax 
any stress above the local slip thresholds, as soon as those slips systems become unstable. 
This relaxation must be finished before we increment the load again. 

To check for unstable slip systems, we can iterate over the elements, and for each element, 
iterate again over its slip systems. If we find an slip with a negative barrier (i.e., the shear 
stress on its slip plane is above its critical value), we add a slip event to its parent element, as
we saw in the previous tutorial. This procedure is fine, but we need to consider what would 
happen if more than a single slip system becomes simultaneously unstable within the same 
element simultaneously. Fortunately, this is no problem here since for the element class 
@ref example::element::Scalar<dim>, only one of its two slip systems can be, by definition, 
unstable. The reason is that they represent the two slip directions of a single slip plane, 
with shear orientations given by the angles \f$ 0 \f$ and \f$ \pi/2 \f$, respectively. 
 
@note the external loading represents a positive and homogeneous contribution to the
xy-component of the stress tensor. However, the local stress is the superposition of the external
stress and the internal stress. The internal stress is induced by the plastic strain field, and is 
heterogeneous with non-zero xx-, yy- and xy- components. Moreover, these components can be negative
in some locations, depending on the spatial position relative to a plastic event. 
Consequently, in general, the magnitude and sign of the local shear stress fluctuates in space 
(and time), activating slips in both the \f$ 0 \f$ and \f$ \pi/2 \f$ orientations.
   
In summary, the process works as follows: we iterate over the elements, and for each element, we 
iterate over its slip systems. We check if the slip's barrier is negative, in which case we add a
slip event to a vector of events and move to the next element. After we have checked all the 
elements, we add the vector of events to the system. The purpose of adding a vector of slip 
events instead of individual events is that the events occur simultaneously. By adding them 
through a vector they will be recorded in the event history in this way. Also, the elastic fields
will be computed only once for all the events simultaneously, which increases the performance a 
great deal. Also, after iterating over the elements, we check if we added any slip event. If we 
didn't, that means that the system is stable everywhere, and we can escape from the relaxation 
loop. If we added at least one new event, the stress fields have been updated. Therefore we must
check again everywhere for instabilities. In this case, we repeat the relaxation loop. This 
process gives rise to a cascade of slip events. The external load increment triggers the first
generation events. The second-generation events are triggered by the first generation ones, and so on.
     
    
```cpp   
      std::vector< mepls::event::Plastic<dim> > events_relax_step;
      bool continue_relaxation = true;

      // relaxation loop (cascade of slip events)
      while( continue_relaxation )
      {
         events_relax_step.clear();

         // whenever we find an unstable slip system, we the corresponding plastic
         // event to the system. Instead of one at a time, we add a vector of them
         // to consider them simoultaneous in time (also, this increases performce since
         // the FEM problem is solved only once for all of them)
         for(auto &element : system)
            for(auto &slip : *element)
               if(slip->barrier < 0.)
               {
                  mepls::event::Plastic<dim> plastic_event(slip);
                  events_relax_step.push_back(plastic_event);

                  // go to the next element, since the other slip system cannot have
                  // a negative barrier too
                  break;                
               }

         system.add(events_relax_step);

         // if some slip events were added, the stress field has changed, and new slip systems
         // might be unstable, so we repeat the relaxation loop to check it. However, if no slip
         // event was added, which means the all the slip systems are stable, and no changes can occur
         // anymore. In this case, don't repeat the relaxation loop
         continue_relaxation = events_relax_step.size() > 0;

      } // relaxation loop
```

### Writing output data {#writing_data}

Now, the system is relaxed, i.e., no more plastic activity takes place until we increase again the 
externally applied strain. At this stable state, we tell the history that the macrostate of the 
system must be recorded. In contrast with the slip and driving events, we are responsible for telling
the history object when the macrostate should be recorded since it does not know 
which state is meaningful to us. In this case, it is only the equilibrium ones.
 
```cpp
      // record the system's macroscale properties. We do it after the relaxation loop
      // so we record a stable state.
      sim_history.add_macro(system);

   } // main simulation loop
   
   // remember that we allocated dynamically using the new operator
   // when we build a model, it will be our responsibility to delete them
   for(auto &element : elements)
      delete element;    
```
   
This two-step dynamics (load increments followed by a relaxation in the form of) will repeat until the 
applied strain reaches the target value of 5%, as established in the condition of the main 
simulation loop. 

When we scape the main simulation loop, we want to save the evolution of some macroscale 
properties. For this, we can iterate over the `macro_evolution` member of the history 
object, @ref mepls::history::History<dim>::macro_evolution. This member is a vector of structs 
@ref mepls::history::History<dim>::MacroSummaryRow, which contain many different macroscale 
properties of the system.

Specifically, we are interested in the applied strain and external stress. Both scalar 
quantities refer to the xy-components of the tensors and will allow us to see plot the 
stress-strain curve of the material. We save the data in CSV format, that is, in two columns 
separated by a comma. The first line contains the name of each column.

```cpp   
   // write in CSV format the total strain and external stress columns from the macroevolution
   // stored in the history object
   std::ofstream output_file("out.csv");

   // columns names
   output_file << "total_strain,ext_stress\n";

   // the data
   for(auto &row : sim_history.macro_evolution)
      output_file << row.total_strain << "," << row.ext_stress << "\n";

   output_file.close();
}   
```
   
   
# Results{#results}

When running the program, the output will consist of two columns. The first one is the total 
applied strain, and the second is the external stress. Instead of showing here the raw output, we 
will redirect it to a file and use a very simple Python script (see next) to plot it:

```sh
$ ./run_sim

$ head out.csv
0.0001,0.006
0.0002,0.012
0.0003,0.018
0.0004,0.024
0.0005,0.03
0.0006,0.036
0.0007,0.042
0.0008,0.048
0.0009,0.054

$ python plot.py out.dat
```

@note when you run `plot.py`, you might need to adapt the command above to use the right path to `out.dat`.

You can find the plotting script `plot.py` in this tutorial's directory. It needs the module 
`pandas`. The script is the following:

```python
import pandas as pd
import matplotlib.pyplot as plt
import sys

# read the data
data = pd.read_csv( sys.argv[1] )

# strain in % and stress in MPa
data['total_strain'] *= 100
data['ext_stress'] *= 1000

data.plot(x='total_strain', y='ext_stress', legend=None)

# tune plot details
plt.xlabel(r'$\varepsilon_{\rm xy}$ (%)', fontsize=15)
plt.ylabel(r'$\Sigma_{\rm xy}$ (MPa)', fontsize=15)
plt.tick_params(labelsize=13)
plt.tight_layout()  
plt.show()
```

The resulting plot is:
<center><img src="step3_stress_strain_curve.png" width="30%"></center>

We see how the stress-strain response adopts the traditional shape observed for most 
materials. Namely, an initial linear regime with a slope defined by the shear modulus, followed 
plastic deformation reflected in the stress drops. The characteristics of this curve depend on the simulation parameters used, and their impact has
been widely studied, see, e.g.,

@note it is common in the literature of elasto-plastic mesoscale models to rescale the stress units 
and work in units of `lambda` or of average slip threshold. The value of the `shear_modulus` is then
given in units of `lambda` instead of GPa. On the other hand, the strains are rescaled by `shear_modulus / lambda`.
This system of units leverages the fact that the system's dynamics
are not sensitive to all the simulation parameters independently but to some adimensional ratio of them.
In this case, that adimensional ratio is `shear_modulus * gamma / lambda`, and it tunes the 
influence of a slip event on its neighborhood (i.e., the intensity of the stress field induced by
the event, of the order of `shear_modulus * gamma`) over the typical resistance `lambda` of the 
neighborhood to imitate it (i.e., to also slip).


In this tutorial, we saw how to control the system according to rules representing a specific 
physical scenario, namely a slowly driven material in athermal conditions. In the next tutorial, we 
will improve the implementation of the same model developed here. To this end, we will MEPLS 
built-in dynamic protocols to implement the same dynamics more efficiently. 

We performed some rudimentary output to obtain the stress-strain curve. This output was very simple 
but was enough for illustrating the tutorial results. In the next tutorial, we will also show how
 to do a full output, which will allow us to visualize the evolution of the system more in detail.
 

# The complete program{#full}

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
#include <mepls/solver.h>
#include <mepls/system.h>
#include <mepls/history.h>
#include <random>

// new header
#include <fstream>

int main()
{
	constexpr unsigned dim = 2;
	std::mt19937 generator(1234567);

	// let's consider a material with a shear modulus of 30 GPa,
	double G = 30.;
	double nu = 0.3;
	dealii::SymmetricTensor<4, dim> C = mepls::utils::tensor::make_isotropic_stiffness<dim>(G, nu);

	unsigned int Nx = 16;
	unsigned int Ny = 16;

	mepls::element::Vector<dim> elements;
	for(double n = 0; n < Nx * Ny; ++n)
	{
		example::element::Scalar<dim>::Config conf;
		conf.number = n;

		// plastic increments induce a local shear deformation of 5%
		conf.gamma = 0.05;

		// the slip thresholds have a scale of 1 GPa (see Weibull dist. for the relation between
		// lambda, k and the average value)
		conf.lambda = 1;

		// the slip thresholds have some disorder. For that, we use a not too-high k (for k->inf
		// the disorder vanishes, while it increases for k->0)
		conf.k = 6;

		auto element = new example::element::Scalar<dim>(conf, generator);
		element->C(C);

		elements.push_back(element);
	}

	mepls::elasticity_solver::LeesEdwards<dim> solver(Nx, Ny);
	for(auto &element : elements)
		solver.set_elastic_properties(element->number(), element->C());
	solver.setup_and_assembly();

	mepls::element::calculate_ext_stress_coefficients(elements, solver);
	mepls::element::calculate_local_stress_coefficients_central(elements, solver);

	mepls::system::Standard<dim> system(elements, solver, generator);


	// we create a history object to record the evolution of the sytem
	mepls::history::History<dim> sim_history("Simulation_history");

	// the system needs to know the history object so that it can record the driving
	// and slip events that we add
	system.set_history(sim_history);


	//----- DYNAMICS OF THE SYSTEM ------

	// main simulation loop till 5% applied shear strain (epsilon_xy, not gamma)
	while(system.macrostate["total_strain"] < 0.05)
	{

		// print some output
		std::cout << system.macrostate["total_strain"] << " " << system.macrostate["ext_stress"]
				  << std::endl;

		//----- PERTURB THE SYTEM ------

		// the applied shear strain is (instantaneously) increased by a small
		// ammout of 0.01%, and is kept fixed afterward
		mepls::event::Driving<dim> driving_event;
		driving_event.dload = 0.0001;
		system.add(driving_event);


		//----- RELAX THE SYTEM ------

		std::vector<mepls::event::Plastic<dim> > events_relax_step;
		bool continue_relaxation = true;

		// relaxation loop (cascade of slip events)
		while(continue_relaxation)
		{
			events_relax_step.clear();

			// whenever we find an unstable slip system, we the corresponding plastic
			// event to the system. Instead of one at a time, we add a vector of them
			// to consider them simoultaneous in time (also, this increases performce since
			// the FEM problem is solved only once for all of them)
			for(auto &element : system)
				for(auto &slip : *element)
					if(slip->barrier < 0.)
					{
						mepls::event::Plastic<dim> plastic_event(slip);
						events_relax_step.push_back(plastic_event);

						// go to the next element, since the other slip system cannot have
						// a negative barrier too
						break;
					}

			system.add(events_relax_step);

			// if some slip events were added, the stress field has changed, and new slip systems
			// might be unstable, so we repeat the relaxation loop to check it. However, if no slip
			// event was added, which means the all the slip systems are stable, and no changes can occur
			// anymore. In this case, don't repeat the relaxation loop
			continue_relaxation = events_relax_step.size() > 0;

		} // relaxation loop

		// record the system's macroscale properties. We do it after the relaxation loop
		// so we record a stable state.
		sim_history.add_macro(system);

	} // main simulation loop

	// remember that we allocated dynamically using the new operator
	// when we build a model, it will be our responsibility to delete them
	for(auto &element : elements)
		delete element;


	// write in CSV format the total strain and external stress columns from the macroevolution
	// stored in the history object
	std::ofstream output_file("out.csv");

	// columns names
	output_file << "total_strain,ext_stress\n";

	// the data
	for(auto &row : sim_history.macro_evolution)
		output_file << row.total_strain << "," << row.ext_stress << "\n";

	output_file.close();
}
```

[deal.II]: https://www.dealii.org/
