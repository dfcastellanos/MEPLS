@page Step4

# Step 4

[< previous ](@ref Step3) | [next >](@ref Step5)

### Table of contents

- [Introducion](#introducion) 
    - [The model](#the_model)
- [The commented program](#comented_program)
    - [Using input parameters files](#param_files)
    - [Writing structured output data](#write_data)
    - [The simulation](#simulation)
    - [The main function](#main)
- [Results](#results)
- [The complete program](#full)


# Introduction{#introducion}

This tutorial will use the same model introduced in the previous tutorial step with little 
modifications. This time, however, we will use MEPLS built-in tools for defining the dynamics of 
the system. Also, we will introduce the last missing parts that the program needs to be user-friendly:
we will show how to set the model parameters using input configuration files, and how to write 
all the simulation the results in a clean and well-structured manner, easy to be read afterward
by Python scripts.


### The model{#the_model}

We consider the same model as in the previous step, i.e., a material being driven in 
strain-controlled pure shear conditions in the athermal quasistatic limit. This time, however, we
consider that the material has some thermal or mechanical history prior to the beginning of the 
driving loading protocol, e.g. the material has undergone aging at non-zero temperature or has 
been strained before.

In this situation, the statistics of the material's structural properties exhibit differences from
the properties renewed just after a local plastic event. To see it more clearly, we can 
think, e.g., of a glass sample annealed during some time. In this situation, the glassy structure 
evolves finding ever more thermodynamically stable configurations. From a micromechanics point of
view, more stable configurations mean higher local slip thresholds. However, when a local 
plastic event takes place, the local atomistic configuration changes, which will result in new 
local slip thresholds and, typically in a structure softer than the initial, highly stable one. 
There are many posibilities to model this kind of behavior, but we will use a fairly simple one. 
To this end, we consider two different slip threshold scale parameters, \f$ \lambda_{\rm init} \f$ 
for 
the intiaal state, and  \f$ \lambda_{\rm renew} \f$ for renewing the thresholds after plastic slip
events occur. Thus, if \f$ \lambda_{\rm init} = \lambda_{\rm renew}\f$ the model reduces entirely
to the previous tutorial step model. However, if \f$ \lambda_{\rm init} > \lambda_{\rm renew}\f$,
 we model a material with an initial structure which is more stable than after plastic 
 deformation. We can 
think of it as a kind of softening, and as we will see in the results section, it leads to a 
stress overshoot in the stress-strain curve and to strain localization in the form of a permanent shear 
band.



# The commented program{#comented_program}
   
We include the same headers of the previous tutorial (see @ref Step3), plus the headers with the tools
for the MEPLS system dynamics, loading input parameter files, and creating output JSON files.
   
```cpp
#include <example.h>
#include <mepls/utils.h>
#include <mepls/solver.h>
#include <mepls/system.h>
#include <mepls/history.h>

// MEPLS built-in dynamics
#include <mepls/dynamics.h>

// to parse command line arguments
#include <cmdparser.hpp>

// to parse input parameters files
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/path_search.h>

// to save output data in JSON format
#include <boost/property_tree/json_parser.hpp>
```

### Using input parameters files{#param_files}

In this tutorial, we want to pass the simulation parameters using a text file. This text file 
will be parsed, and the values will be set in a struct of type `Parameters`  that we will use
throughout the program. The `Parameters` struct uses 
[deal.II's ParameterHandler](https://dealii.org/current/doxygen/deal.II/classParameterHandler.html) class 
for parsing the input text file. In the definition of the struct, we first declare the different 
members and their default values. Then, in `declare_entries(...)`, we 
declare the entries that the parser will look for in the input file. The definition of 
each entry needs a default value (the second argument of the function), for which we reuse the 
default values defined before. After that, in `load_entries(...)` we define
how each member gets its value from the entries of the parsed file.

```cpp
struct Parameters
{
   // these are the parameters and we will use in the simulation, initialized with some default values
   unsigned int seed = 1234567;
   unsigned int Nx = 32;
   unsigned int Ny = 32;
   double G = 30.;
   double nu = 0.3;
   double gamma = 0.05;
   double k = 6.;
   double strain_limit = 0.05;
   double lambda_init = 1.;
   double lambda_renew = 1.;
   std::string filename = "out.json";

   void declare_entries(dealii::ParameterHandler &prm)
   {
      // We declare the entries of the parameters text file. Each entry matches the name of a
      // simulation parameters (although it doesn't need to) and has a default value. We use as
      // the default value is the same value of the variables declared above

      prm.enter_subsection("Section1");

      prm.declare_entry("seed", mepls::utils::str::to_string(seed), dealii::Patterns::Integer(0), "");
      prm.declare_entry("Nx", mepls::utils::str::to_string(Nx), dealii::Patterns::Integer(0), "");
      prm.declare_entry("Ny", mepls::utils::str::to_string(Ny), dealii::Patterns::Integer(0), "");
      prm.declare_entry("G", mepls::utils::str::to_string(G), dealii::Patterns::Double(0.0), "");
      prm.declare_entry("nu", mepls::utils::str::to_string(nu), dealii::Patterns::Double(0.0), "");
      prm.declare_entry("gamma", mepls::utils::str::to_string(gamma), dealii::Patterns::Double(0.0), "");
      prm.declare_entry("strain_limit", mepls::utils::str::to_string(strain_limit), dealii::Patterns::Double(0.0), "");
      prm.declare_entry("lambda_init", mepls::utils::str::to_string(lambda_init), dealii::Patterns::Double(0.0), "");
      prm.declare_entry("lambda_renew", mepls::utils::str::to_string(lambda_renew), dealii::Patterns::Double(0.0), "");
      prm.declare_entry("k", mepls::utils::str::to_string(k), dealii::Patterns::Double(0.0), "");
      prm.declare_entry("filename", filename, dealii::Patterns::FileName(), "");

      prm.leave_subsection();
   }

   void load_entries(dealii::ParameterHandler &prm)
   {
      prm.enter_subsection("Section1");

      // We define how each parameter gets its value from a parameters file entry

      seed = prm.get_integer("seed");
      Nx = prm.get_integer("Nx");
      Ny = prm.get_integer("Ny");
      G = prm.get_double("G");
      nu = prm.get_double("nu");
      gamma = prm.get_double("gamma");
      k = prm.get_double("k");
      strain_limit = prm.get_double("strain_limit");
      lambda_init = prm.get_double("lambda_init");
      lambda_renew = prm.get_double("lambda_renew");
      filename = prm.get("filename");

      prm.leave_subsection();
   }
```

Now, in `load_file(...)` we define how the `Parameters` struct loads the 
input text file. After that, in `generate_file(...)`, we define how to 
generate a new file, that can be used as a template input file.  
  
```cpp
   void load_file(const std::string &filename)
   {
      // load a parameters file

      dealii::ParameterHandler prm;
      declare_entries(prm);
      prm.parse_input(filename);
      load_entries(prm);
   }

   void generate_file(const std::string &filename)
   {
      // generate a parameters file template

      std::ofstream outfile(filename);
      dealii::ParameterHandler prm;
      declare_entries(prm);
      prm.print_parameters(outfile, dealii::ParameterHandler::OutputStyle::Text);
   }
   
}; // Parameters
```

### Writing structured output data{#write_data}
  
we will save the results into a tree data structure of class `boost::property_tree::ptree`. This 
tree allows us to organize the data in a hierarchical way easily. The hierarchy can then be 
easily written to a JSON plain text file. First, we write some metadata, such as a simulation or 
model name and a description. Here we can add any fields we may deem necessary. Then, we save the
simulation parameters. Although they are already stored in the input parameters file, it is 
advisable to save them together with the data, so the final output file is as self-descriptive as
 possible.

Which data to save depends on what we are interested in and the storage capacity. Here, we will 
use the plastic and driving event histories to create the stress-strain curve. These datasets 
will allow us to build the stress-strain curve in combination with spatial 2D maps of plastic 
activity.
 
The tree structure is very flexible, and we can store the data in many different ways. An 
easy and efficient way is to consider the tree as a file system, through which we can navigate to
locate the datasets, which contain the actual heavy data formatted in CSV. We will 
store the datasets as `Data/plastic_events` and `Data/driving_events`. (We could make the 
tree even more descriptive and have metadata associated with individual datasets by doing, e.g.,
`Data/driving_events/metadata` and `Data/driving_events/csv`).  
  
```cpp
template<int dim>
void write_data(
	const mepls::history::History<dim> &sim_history, const Parameters &p)
{
	boost::property_tree::ptree data_tree;

	// write some metadata, such as a simulation name and description
	data_tree.put("Name", "Step4");
	data_tree.put("Description", "System driven in the athermal quasistatic limit");

	// write the simulation parameters
	data_tree.put("Parameters.dim", 2);
	data_tree.put("Parameters.seed", p.seed);
	data_tree.put("Parameters.Nx", p.Nx);
	data_tree.put("Parameters.Ny", p.Ny);
	data_tree.put("Parameters.G", p.G);
	data_tree.put("Parameters.nu", p.nu);
	data_tree.put("Parameters.gamma", p.gamma);
	data_tree.put("Parameters.lambda_renew", p.lambda_renew);
	data_tree.put("Parameters.lambda_init", p.lambda_init);
	data_tree.put("Parameters.k", p.k);
	data_tree.put("Parameters.strain_limit", p.strain_limit);


	// write some metadata, such as a simulation name and description
	data_tree.put("Name", "Step4");
	data_tree.put("Description", "System driven in the athermal quasistatic limit");

	// write the simulation parameters
	data_tree.put("Parameters.dim", 2);
	data_tree.put("Parameters.seed", p.seed);
	data_tree.put("Parameters.Nx", p.Nx);
	data_tree.put("Parameters.Ny", p.Ny);
	data_tree.put("Parameters.G", p.G);
	data_tree.put("Parameters.nu", p.nu);
	data_tree.put("Parameters.gamma", p.gamma);
	data_tree.put("Parameters.lambda_renew", p.lambda_renew);
	data_tree.put("Parameters.lambda_init", p.lambda_init);
	data_tree.put("Parameters.k", p.k);
	data_tree.put("Parameters.strain_limit", p.strain_limit);

	// We write the event histories with CSV format to a string.
	// Here, we write only the columns of interest, but there are more available (see the
	// documentation).
	std::ostringstream plastic_events_csv;
	plastic_events_csv << "index,element,eigenstrain_00,eigenstrain_11,eigenstrain_01\n";
	for(auto &row : sim_history.plastic)
		plastic_events_csv << row.index << "," << row.element << "," << row.eigenstrain_00 << ","
						   << row.eigenstrain_11 << "," << row.eigenstrain_01 << "\n";

	data_tree.put("Data.plastic_events", plastic_events_csv.str());


	std::ostringstream driving_events_csv;
	driving_events_csv << "index,dext_stress,dtotal_strain\n";
	for(auto &row : sim_history.driving)
		driving_events_csv << row.index << "," << row.dext_stress << "," << row.dtotal_strain
						   << "\n";

	data_tree.put("Data.driving_events", driving_events_csv.str());


	std::ostringstream macro_evolution_csv;
	macro_evolution_csv << "index,ext_stress,total_strain,time,av_vm_stress,av_vm_plastic_strain\n";
	for(auto &row : sim_history.macro_evolution)
		macro_evolution_csv << row.index << "," << row.ext_stress << "," << row.total_strain << ","
							<< row.time << "," << row.av_vm_stress << ","
							<< row.av_vm_plastic_strain << "\n";

	data_tree.put("Data.macro_evolution", macro_evolution_csv.str());


	// write the data tree to a JSON file
	std::ofstream output_file(p.filename);
	boost::property_tree::json_parser::write_json(output_file, data_tree);
	output_file.close();
}
``` 

### The simulation{#simulation}

We implement a `run()` function that constitutes the simulation itself, and will be called later 
from the `%main()` function. The `run()` function takes as an input the parameters struct and 
a structure to store the simulation results that we will describe later. In the previous 
tutorials, we have already described in detail the parts of this function. First, we create 
and set up the elements, solver, and system. This time, the values of the parameters are set using
the input parameters struct. Also, we store the elements in two different vectors. The vector 
of class @ref mepls::element::Vector<dim> stores pointers to the element base class @ref 
mepls::element::Element<dim>, which is enough for passing it as an argument to the different MEPLS
tools. However, in this tutorial, we will reaccess the elements' internal configuration struct. To
 do this, we need pointers to the derive class @ref example::element::Scalar<dim>. There are two 
possibilities: we can cast the pointers to the known derive class, or we can store the elements 
also in a vector `std::vector<example::element::Scalar<dim> *>`. We will use the second 
option. Note that havint two vector of elemnts doesn't have any real impact on the 
simulation since the vectors contain only pointers to the actual element objects.

```cpp
void run(const Parameters &p, boost::property_tree::ptree &data_tree)
{
   // we do the same as in the previous tutorial, but this time we use the
   // parameters from the input Parameters object
  
   constexpr unsigned dim = 2;
   std::mt19937 generator(p.seed);

   dealii::SymmetricTensor<4, dim> C = mepls::utils::tensor::make_isotropic_stiffness<dim>(p.G, p.nu);

   mepls::element::Vector<dim> elements;

   // we also store the elements into a vector that knows the derived class, so we can access
   // the example::element::Scalar<dim>::conf struct later (this member cannot be accessed
   // through the pointers to the base mepls::element::Element<dim>)
   std::vector<example::element::Scalar<dim> *> elements_scalar;

   for(double n = 0; n < p.Nx * p.Ny; ++n)
   {
      example::element::Scalar<dim>::Config conf;
      conf.number = n;
      conf.gamma = p.gamma;
      conf.lambda = p.lambda_init;
      conf.k = p.k;

      auto element = new example::element::Scalar<dim>(conf, generator);
      element->C(C);

      elements_scalar.push_back(element);
      elements.push_back(element);
   }

   mepls::elasticity_solver::LeesEdwards<dim> solver(p.Nx, p.Ny);
   for(auto &element : elements)
      solver.set_elastic_properties(element->number(), element->C());
   solver.setup_and_assembly();

   mepls::element::calculate_ext_stress_coefficients(elements, solver);
   mepls::element::calculate_local_stress_coefficients_central(elements, solver);

   mepls::system::Standard<dim> system(elements, solver, generator);

   mepls::history::History<dim> sim_history("Simulation_history");

   system.set_history(sim_history);
```

When the elements are created, their slip systems are initialized using the parameters passed in 
the struct @ref example::element::Scalar<dim>::Config. As explained in the introduction, we want 
that after local plastic deformation has occurred, the slip systems have different statistical 
properties. As seen in @ref Step1, the elements will keep their initial slip systems until their 
structural properties are first renewed. A way to do this is to change the values of the elements 
configuration parameters after they have been created, but before the first structural renewal 
takes place (which occurs in the dynamics, when plastic events are added). To this end, we 
iterate over the elements using the pointers to the derived class, access their configuration 
parameters, and replace the value of the slip threshold scale from the initial `lambda_init` to
 `lambda_renew`.

```cpp
   // when the elements were created, their slip systems were initialized with thresholds from a
   // Weibull distribution with scale lambda_init. We change it now to lambda_renew. Therefore,
   // when plastic deformation occurs, the renew local thresholds will have a different
   // average 
   for(auto &element : elements_scalar)
      element->conf.lambda = p.lambda_renew;
```

Now we implement the dynamical evolution of the system. We use the same rules described 
in the previous tutorial @ref Step3. However, now we use MEPLS built-in functions, which offer an
optimized implementation, debug flags, and extra control through the function arguments.

To check if the simulation must stop, we use the class @ref mepls::utils::ContinueSimulation. When we use 
the call operator, the object will tell us whether the simulation should continue running or not.  
This class can be used to check for conditions in different parts of the simulation, and we can 
associate a different message to each of these conditions. Whenever one the conditions is not 
met, the object's internal state is set to `false`, and it won't change again. At the end, 
using the associated message, the user will be informed about which of the conditions caused 
the simulation to stop. In this case, we will use as the stopping criterion reaching a maximum total 
strain value. Also, passing the continue_simulation object to the function @ref 
mepls::dynamics::relaxation, we will also check for an extra condition within the function body.
Specifically, it will check if the avalanche size overcomes a certain maximum upper limit, in which 
case two possibilities exist: the avalanche is unrealistically big because we are using simulation 
parameters that are far from the reality, or (in some models) such avalanche can mean that a 
material has undergone mechanical failure. When the main loop finishes, we write the data to the 
output file.
   
```cpp
   // This object will allow us to check for different conditions by which
   // the simulation might stop. Different conditions might be checked, but as long as one of 
   // them evaluates to false, when continue_simulation() is called it will return false
   mepls::utils::ContinueSimulation continue_simulation;

   sim_history.add_macro(system);

   // run while the continue_simulation object says so
   while(continue_simulation())
   {
      std::cout << system.macrostate["total_strain"] << " " << system.macrostate["ext_stress"] << std::endl;


      // these are the same dynamics we implemented in the previous tutorial, but using MEPLS built-in

      // apply an external strain increment of 0.01%
      mepls::dynamics::finite_extremal_dynamics_step(1e-4, system);
      sim_history.add_macro(system);

      // perform and avalanche of slip events. By passing the continue_simulation object, 
      // the relaxation function can set its own condition for stopping the simulation. 
      // Specifically, it will check if the avalanche size overcomes a certain maximum upper limit
      mepls::dynamics::relaxation(system, continue_simulation);
      sim_history.add_macro(system);


      // check if the strain has reached the strain limit. If it has, the next time continue_simulation() is called
      // it will return false, so we will exit the main loop
      continue_simulation(system.macrostate["total_strain"] < p.strain_limit, "total strain limit reached");

   }

   // print the message of the stopping condition that was met 
   std::cout << continue_simulation << std::endl;

   for(auto &element : elements)
      delete element;
      
      
	write_data(sim_history, p);
	
} // run
```

### The main function{#main}
   
Finally, we define the `%main()` function. The first thing that we do here is to read the command
line arguments passed to the program. We only define the argument `-f`, which takes the path to 
the input parameters file to be parse. If no `-f` argument is given, a default value of `
./default.prm` will be used.
 
```cpp
int main(int argc, char *argv[])
{
   // Read the command line arguments. We define the -f 'filename' to pass the
   // path to the parameters file
   cli::Parser parser(argc, argv);
   parser.set_optional<std::string>("f", "file", "./default.prm", "Name of the input configuration file");
   parser.run_and_exit_if_error();
   // you can check https://github.com/FlorianRappl/CmdParser for a further documentation of cli::Parser
```

Now we create a `Parameters` struct, introduced at the beginning of this section. We use a `try` 
block to catch the exception thrown if the input file does not exist. If this is the case, a new 
parameters file with that name is generated, and the program ends.

```cpp
   // Create the parameters object
   Parameters p;
   
   // We try to load the parameters file, but if it doesn't exist, we generate a new one with the name
   // default.prm and default values
   try
   {
      p.load_file(parser.get<std::string>("f"));
   }
   catch(dealii::PathSearch::ExcFileNotFound &)
   {
      p.generate_file(parser.get<std::string>("f"));
      std::cout << "Configuration file " << parser.get<std::string>("f") << " created" << std::endl;
      return 1;
   }
```

Finally, we call the `run()` function with the parameters struct as the argument. When it 
returns, the program will end and the simulation results will be stored in a JSON text file. In 
the Results section will show how to load and visualize the data stored in that file. 

```cpp

	// run the simulation
   run(p);

   return 0;

} // main
```


# Results{#results}

After compiling the program (see @ref HowToBuild), we can run it as `./run_sim`,

```sh
$ ./run_sim
Configuration file ./default.prm created
```

Since no parameters file was given and `default.prm` was not found, the program created a 
template file named `default.prm`. We can edit it with a text editor and set the values of the 
parameters we will use during the simulation run. In this example, we will use a parameters 
file like this,

```
# Listing of Parameters
# ---------------------
subsection Section1
  set G             = 30.00
  set Nx            = 32
  set Ny            = 32
  set filename      = out.json
  set gamma         = 0.05
  set k             = 6.00
  set lambda_init   = 1.00
  set lambda_renew  = 1.00
  set nu            = 0.30
  set seed          = 1234567
  set strain_limit  = 0.05
end
```

Now, if we run it again, the file `default.prm` will be used. We could also pass 
a different file as `./run_sim -f other_file.prm`. When the simulation ends, the program will 
create a file `out.json` containing the results. A basic way of loading the data generated by
this program using Python might look as follows,

```python
import json
from io import StringIO
import pandas as pd

sim_data = json.load( open('out.json') )

df = pd.read_csv( StringIO(sim_data['Data']['plastic_events']) )
```

Based on this approach for loading the data, the Python script `animate.py` located in 
the `MEPLS/python` directory creates an interactive animation of the results. You can run
the script as follows:

```bash
python animate.py --data out.json /Data/plastic_events /Data/driving_events --vars total_strain ext_stress dstrain --rescale 100 1000 1 --labels '$\varepsilon_{\rm xy}$ (%)' '$\Sigma_{\rm xy}$ (MPa)' '$\varepsilon_{\rm vm}(\vec{r})$' 
```

The first argument indicates that we are reading the file `out.json` (you need to use the right 
path to the file `out.json`). The second and third that from that file we load the plastic and 
driving event histories located in `/Data/plastic_events` and `/Data/driving_events` respectively
(this paths are the same as defined by us in the program when filling in the data tree). The 
animation consists of two plots, one with a curve Y(x) and another with a 2D plot of a variable Z. 
The argument forth, fifth and sixth arguments indicate, respsectively, the names of X, Y and Z. 
With the command above, we are going to plot the `ext_stress` vs. `total_strain` curve along with 
a 2D plot of `dstrain` (which, in the script, is a shortcut for the local von Mises plastic
strain increments). The rest of the arguments are optional. The argument `--rescale` means that we
 are going rescale the `total_strain` by 100, the `ext_stress` by 1000 and the `dstrain` remains 
 unchanged. This rescaling operation is transforming the stress from GPa to MPa and the strain 
 into a percentage. The argument `--labels` means that we are going to use the axis labels \f$ 
\varepsilon_{\rm xy} \rm{(\%)} \f$, \f$ \Sigma_{\rm xy} \rm{(MPa)} \f$ and \f$ \varepsilon_{\rm vm}(\vec{r})\f$. 
The resulting animation has interactive controls, which you can find with `python animate.py 
--help`. 

The following GIFs, showing the evolution of the material's stress-strain response and the spatial map of 
plastic activity (specifically, the local von Mises plastic strain) are created from the animations. 

For the simulation performed with `lambda_init = lambda_renew = 1.0`:
<center><img src="step4_stress_strain_curve_flat.gif" width="40%"></center>

This the simulation performed with with `lambda_init = 1.3` and `lambda_renew = 1.0`:
<center><img src="step4_stress_strain_curve_drop.gif" width="40%"></center>

In the first case (`lambda_init = lambda_renew = 1.0`), the results are the same as in @ref 
Step3. Initially, we find the elastic regime after which the system reaches the stationary flow 
regime, where the external stress remains on average constant. The spatial map initially shows 
sparse activity since the activation of slip events is dominated by the least stable slip 
systems, whose locations are random and (initially) uncorrelated. However, after some plastic 
deformation, an internal stress field build-up, which induces spatial correlations. The effect of
correlations can be observed in how the deformation localizes in transient horizontal and 
vertical bands. Note that, despite this localization, the deformation is still statistically 
homogeneous since the loading conditions induce a spatially-homogeneous shear stress field. 

However, in the second case (`lambda_init = 1.3` and `lambda_renew = 1.0`), the slip
thresholds' scale is higher in the initial state than in the renewed one (i.e., after local plastic 
deformation). We can observe how this has profound implications in the 
response of the material. Specifically, we observe a stress overshoot followed by a significant stress 
drop. The spatial map shows the formation of a macroscopic permanent shear band with little 
plastic deformation outside of it. For brittle materials, such a shear band can lead to rupture. The 
stress overshoot and how abrupt the stress drop is will depend on the difference 
between of `lambda_init` and `lambda_renew`. The other simulation parameters also have an 
impact on it as, e.g., the parameter `k` controlling the disorder in the slip thresholds. The 
effects of the different parameters on the stress overshoot and the formation of the shear band 
have been widely studied in the scientific literature, see e.g., .

In this tutorial, we implemented the system dynamics using MEPLS built-in tools. Also, we saw how
to control the simulation parameters using input files and how to create and load output files 
with complex structured data on them. At this point, we have described all the main tools 
necessary to develop a model and simulate it. 

In the following tutorial, we will consider an initial state with properties that differ from the
renewed ones, as in this tutorial. However, the initial state will be dynamically constructed
by simulating a sample deforming under creep conditions at some non-zero temperature. Also, we 
will see how run several instances of a simulation in parallel, which will improve the 
performance of the statistical sampling of the model.



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

// MEPLS built-in dynamics
#include <mepls/dynamics.h>

// to parse command line arguments
#include <cmdparser.hpp>

// to parse input parameters files
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/path_search.h>

// to save output data in JSON format
#include <boost/property_tree/json_parser.hpp>


struct Parameters
{
	// these are the parameters and we will use in the simulation, initialized with some default values
	unsigned int seed = 1234567;
	unsigned int Nx = 32;
	unsigned int Ny = 32;
	double G = 30.;
	double nu = 0.3;
	double gamma = 0.05;
	double k = 6.;
	double strain_limit = 0.05;
	double lambda_init = 1.;
	double lambda_renew = 1.;
	std::string filename = "out.json";

	void declare_entries(dealii::ParameterHandler &prm)
	{
		// We declare the entries of the parameters text file. Each entry matches the name of a
		// simulation parameters (although it doesn't need to) and has a default value. We use as
		// the default value is the same value of the variables declared above

		prm.enter_subsection("Section1");

		prm.declare_entry("seed", mepls::utils::str::to_string(seed), dealii::Patterns::Integer(0),
						  "");
		prm.declare_entry("Nx", mepls::utils::str::to_string(Nx), dealii::Patterns::Integer(0), "");
		prm.declare_entry("Ny", mepls::utils::str::to_string(Ny), dealii::Patterns::Integer(0), "");
		prm.declare_entry("G", mepls::utils::str::to_string(G), dealii::Patterns::Double(0.0), "");
		prm.declare_entry("nu", mepls::utils::str::to_string(nu), dealii::Patterns::Double(0.0), "");
		prm.declare_entry("gamma", mepls::utils::str::to_string(gamma), dealii::Patterns::Double(0.0), "");
		prm.declare_entry("strain_limit", mepls::utils::str::to_string(strain_limit), dealii::Patterns::Double(0.0), "");
		prm.declare_entry("lambda_init", mepls::utils::str::to_string(lambda_init), dealii::Patterns::Double(0.0), "");
		prm.declare_entry("lambda_renew", mepls::utils::str::to_string(lambda_renew), dealii::Patterns::Double(0.0), "");
		prm.declare_entry("k", mepls::utils::str::to_string(k), dealii::Patterns::Double(0.0), "");
		prm.declare_entry("filename", filename, dealii::Patterns::FileName(), "");

		prm.leave_subsection();
	}

	void load_entries(dealii::ParameterHandler &prm)
	{
		prm.enter_subsection("Section1");

		// We define how each parameter gets its value from a parameters file entry

		seed = prm.get_integer("seed");
		Nx = prm.get_integer("Nx");
		Ny = prm.get_integer("Ny");
		G = prm.get_double("G");
		nu = prm.get_double("nu");
		gamma = prm.get_double("gamma");
		k = prm.get_double("k");
		strain_limit = prm.get_double("strain_limit");
		lambda_init = prm.get_double("lambda_init");
		lambda_renew = prm.get_double("lambda_renew");
		filename = prm.get("filename");

		prm.leave_subsection();
	}

	void load_file(const std::string &filename)
	{
		// load a parameters file

		dealii::ParameterHandler prm;
		declare_entries(prm);
		prm.parse_input(filename);
		load_entries(prm);
	}

	void generate_file(const std::string &filename)
	{
		// generate a parameters file template

		std::ofstream outfile(filename);
		dealii::ParameterHandler prm;
		declare_entries(prm);
		prm.print_parameters(outfile, dealii::ParameterHandler::OutputStyle::Text);
	}

};


template<int dim>
void write_data(
	const mepls::history::History<dim> &sim_history, const Parameters &p)
{
	boost::property_tree::ptree data_tree;

	// write some metadata, such as a simulation name and description
	data_tree.put("Name", "Step4");
	data_tree.put("Description", "System driven in the athermal quasistatic limit");

	// write the simulation parameters
	data_tree.put("Parameters.dim", 2);
	data_tree.put("Parameters.seed", p.seed);
	data_tree.put("Parameters.Nx", p.Nx);
	data_tree.put("Parameters.Ny", p.Ny);
	data_tree.put("Parameters.G", p.G);
	data_tree.put("Parameters.nu", p.nu);
	data_tree.put("Parameters.gamma", p.gamma);
	data_tree.put("Parameters.lambda_renew", p.lambda_renew);
	data_tree.put("Parameters.lambda_init", p.lambda_init);
	data_tree.put("Parameters.k", p.k);
	data_tree.put("Parameters.strain_limit", p.strain_limit);


	// write some metadata, such as a simulation name and description
	data_tree.put("Name", "Step4");
	data_tree.put("Description", "System driven in the athermal quasistatic limit");

	// write the simulation parameters
	data_tree.put("Parameters.dim", 2);
	data_tree.put("Parameters.seed", p.seed);
	data_tree.put("Parameters.Nx", p.Nx);
	data_tree.put("Parameters.Ny", p.Ny);
	data_tree.put("Parameters.G", p.G);
	data_tree.put("Parameters.nu", p.nu);
	data_tree.put("Parameters.gamma", p.gamma);
	data_tree.put("Parameters.lambda_renew", p.lambda_renew);
	data_tree.put("Parameters.lambda_init", p.lambda_init);
	data_tree.put("Parameters.k", p.k);
	data_tree.put("Parameters.strain_limit", p.strain_limit);

	// We write the event histories with CSV format to a string.
	// Here, we write only the columns of interest, but there are more available (see the
	// documentation).
	std::ostringstream plastic_events_csv;
	plastic_events_csv << "index,element,eigenstrain_00,eigenstrain_11,eigenstrain_01\n";
	for(auto &row : sim_history.plastic)
		plastic_events_csv << row.index << "," << row.element << "," << row.eigenstrain_00 << ","
						   << row.eigenstrain_11 << "," << row.eigenstrain_01 << "\n";

	data_tree.put("Data.plastic_events", plastic_events_csv.str());


	std::ostringstream driving_events_csv;
	driving_events_csv << "index,dext_stress,dtotal_strain\n";
	for(auto &row : sim_history.driving)
		driving_events_csv << row.index << "," << row.dext_stress << "," << row.dtotal_strain
						   << "\n";

	data_tree.put("Data.driving_events", driving_events_csv.str());


	std::ostringstream macro_evolution_csv;
	macro_evolution_csv << "index,ext_stress,total_strain,time,av_vm_stress,av_vm_plastic_strain\n";
	for(auto &row : sim_history.macro_evolution)
		macro_evolution_csv << row.index << "," << row.ext_stress << "," << row.total_strain << ","
							<< row.time << "," << row.av_vm_stress << ","
							<< row.av_vm_plastic_strain << "\n";

	data_tree.put("Data.macro_evolution", macro_evolution_csv.str());


	// write the data tree to a JSON file
	std::ofstream output_file(p.filename);
	boost::property_tree::json_parser::write_json(output_file, data_tree);
	output_file.close();
}


void run(const Parameters &p)
{
	//----- SETUP ------

	// we do the same as in the previous tutorial, but this time we use the
	// parameters from the input Parameters object

	constexpr unsigned dim = 2;
	std::mt19937 generator(p.seed);

	dealii::SymmetricTensor<4, dim> C = mepls::utils::tensor::make_isotropic_stiffness<dim>(p.G,
																							p.nu);

	mepls::element::Vector<dim> elements;

	// we also store the elements into a vector that knows the derived class, so we can access
	// the example::element::Scalar<dim>::conf struct later (this member cannot be accessed
	// through the pointers to the base mepls::element::Element<dim>)
	std::vector<example::element::Scalar<dim> *> elements_scalar;

	for(double n = 0; n < p.Nx * p.Ny; ++n)
	{
		example::element::Scalar<dim>::Config conf;
		conf.number = n;
		conf.gamma = p.gamma;
		conf.lambda = p.lambda_init;
		conf.k = p.k;

		auto element = new example::element::Scalar<dim>(conf, generator);
		element->C(C);

		elements_scalar.push_back(element);
		elements.push_back(element);
	}

	mepls::elasticity_solver::LeesEdwards<dim> solver(p.Nx, p.Ny);
	for(auto &element : elements)
		solver.set_elastic_properties(element->number(), element->C());
	solver.setup_and_assembly();

	mepls::element::calculate_ext_stress_coefficients(elements, solver);
	mepls::element::calculate_local_stress_coefficients_central(elements, solver);

	mepls::system::Standard<dim> system(elements, solver, generator);

	mepls::history::History<dim> sim_history("Simulation_history");

	system.set_history(sim_history);

	// when the elements were created, their slip systems were initialized with thresholds from a
	// Weibull distribution with scale lambda_init. We change it now to lambda_renew. Therefore,
	// when plastic deformation occurs, the renew local thresholds will have a different
	// average
	for(auto &element : elements_scalar)
		element->conf.lambda = p.lambda_renew;


	//----- DYNAMICS ------

	// This object will allow us to check for different conditions by which
	// the simulation might stop. Different conditions might be checked, but as long as one of
	// them evaluates to false, when continue_simulation() is called it will return false
	mepls::utils::ContinueSimulation continue_simulation;

	sim_history.add_macro(system);

	// run while the continue_simulation object says so
	while(continue_simulation())
	{
		std::cout << system.macrostate["total_strain"] << " " << system.macrostate["ext_stress"]
				  << std::endl;


		// these are the same dynamics we implemented in the previous tutorial, but using MEPLS built-in

		// apply an external strain increment of 0.01%
		mepls::dynamics::finite_extremal_dynamics_step(1e-4, system);
		sim_history.add_macro(system);

		// perform and avalanche of slip events. By passing the continue_simulation object,
		// the relaxation function can set its own condition for stopping the simulation.
		// Specifically, it will check if the avalanche size overcomes a certain maximum upper limit
		mepls::dynamics::relaxation(system, continue_simulation);
		sim_history.add_macro(system);


		// check if the strain has reached the strain limit. If it has, the next time continue_simulation() is called
		// it will return false, so we will exit the main loop
		continue_simulation(system.macrostate["total_strain"] < p.strain_limit,
							"total strain limit reached");

	}

	// print the message of the stopping condition that was met
	std::cout << continue_simulation << std::endl;

	for(auto &element : elements)
		delete element;


	write_data(sim_history, p);
}


int main(int argc, char *argv[])
{
	// Read the command line arguments. We define the -f 'filename' to pass the
	// path to the parameters file
	cli::Parser parser(argc, argv);
	parser.set_optional<std::string>("f", "file", "./default.prm",
									 "Name of the input configuration file");
	parser.run_and_exit_if_error();
	// you can check https://github.com/FlorianRappl/CmdParser for a further documentation of cli::Parser



	// Create the parameters object
	Parameters p;

	// We try to load the parameters file, but if it doesn't exist, we generate a new one with the name
	// default.prm and default values
	try
	{
		p.load_file(parser.get<std::string>("f"));
	}
	catch(dealii::PathSearch::ExcFileNotFound &)
	{
		p.generate_file(parser.get<std::string>("f"));
		std::cout << "Configuration file " << parser.get<std::string>("f") << " created"
				  << std::endl;
		return 1;
	}


	// run the simulation
	run(p);


	return 0;
}
```

[deal.II]: https://www.dealii.org/