@page Step5

# Step 5

[< previous ](@ref Step4) | >

### Table of contents

- [Introducion](#introducion_5) 
    - [The model](#the_model_5)
    - [The Kinetic Monte Carlo method for thermal slip activation](#kmc)
    - [Setting the temperature](#setting_temperature)
- [The commented program](#comented_program)
    - [Creep simulation](#creep)
    - [Transition step and AQS simulation](#transition)
- [Results](#results_5)
- [The complete program](#full_5)


# Introduction{#introducion_5}

This tutorial will expand the model considered in the two previous tutorials. We will see how 
different driving protocols can be connected, which allows us to simulate a material's history in a
complete manner. Specifically, we will run an AQS simulation using a 
system's initial configuration that results from a previous simulation under creep conditions. Also, we will use 
thread-based parallelism to run different repetitions of the same simulation simultaneously.

## The model{#the_model_5}

This tutorial will continue with the model of @ref Step3 and @ref Step4, namely a material 
sample driving in the athermal quasistatic limit. In @ref Step4, we considered that the system's 
initial configuration had slip thresholds, on average, higher than once they are renewed. We did 
that by explicitly considering two different slip threshold scale parameters. We saw that this 
choice lead to strain localization into a shear band and a stress overshoot in the stress-strain 
curve. 

Here, we want to reproduce the same kind of behavior in a more physically grounded manner. To this 
end, instead of explicitly considering different thresholds scale parameters, we will have only 
one parameter. However, before the AQS simulation, we will simulate a sample deforming under 
creep conditions, i.e., at non-zero temperature and constant external stress 
@cite Castellanos2018 @cite Castellanos2019. During the creep process, the slip thresholds will 
undergo a survival bias, by which lower thresholds tend to be renewed by higher ones (see the results 
section of @ref Step3 for a short discussion about this). Due to such bias, after some time, the
existing thresholds will be on average higher than as-renewed from their Weibull pdf. 
Consequently, when using such a state as the initial one in the AQS simulation, the AQS 
stress-strain curve will exhibit a behavior similar to @ref Step4. However, in contrast with @ref
 Step4, the magnitude of the stress overshoot and the intensity of the strain localization will 
 be a direct result of a well-defined sample's thermal and mechanical history. 


## The Kinetic Monte Carlo method for thermal slip activation{#kmc}

To simulate thermal activation of plastic activity during the creep simulation, we use the Kinetic 
Monte Carlo (KMC) method @cite DFCastellanos_CRP @cite Castellanos2019 @cite Castellanos2018 
@cite FernandezCastellanos2019. This method requires the knowledge of all the possible transitions 
 that the system can make from the current state towards a new one, and the energy barrier 
 associated with each transition.

In this case, each possible transition corresponds to the activation of a specific slip system. 
The energy barrier \f$ \Delta E \f$ for a specific slip activation can be related to its stress 
distance to threshold (@ref mepls::slip::Slip<dim>::barrier) \f$ \Delta \tau^{\rm c} \f$ as \f$ 
\Delta E \approx \Delta \tau^{\rm c} V_{\rm a}\f$. The quantity \f$ V_{\rm a} \f$ is the so-called 
activation volume, which is of the order of the product of the typical local strain induced by a 
plastic event and the the volume occupied by the event. It is a microscopic quantity 
characteristic of a specific microstructure, and is an input to the model.

The KMC method models thermal activation as a Poisson process, where each possible transition (i
.e., activation of a slip system) is an independent Poisson variable with an activation rate \f$ 
\nu \f$. Here, we consider an Arrenius dependency on temperature and energy,

\f[ \nu(n) = \nu_0 e^{-\frac{\Delta E(n)}{k_{\rm B}T}}\f]

where \f$ n \f$ denotes a specific slip sytem and the parameter \f$ \nu_0 \f$ is a microscopic slip 
activation rate characteristic of the microstructure. The units of time in the model are defined
 by \f$ \Delta t_0 = 1/\nu_0\f$. In terms of stress, we can approximate the expression above as,

\f[ \nu(n) = \nu_0 e^{-\frac{\Delta\tau^{\rm c}(n) V_{\rm a}}{k_{\rm B}T}} = \nu_0 e^{-\frac{\Delta\tau^{\rm c}(n)}{T^{\prime}}} \f]

where \f$ T^{\prime} = k_{\rm B}T/V_{\rm a} \f$. Note that this rescaled temperature has units of
stress, and it characterizes the typical amplitude of the local stress fluctuations induced by 
temperature. 

In summary, we assign to each slip system \f$ n \f$ an activation rate \f$ \nu(n) \f$ based on 
its distance to threshold \f$ \Delta \tau^{\rm c}(n) \f$. Then we simulate the sequence of 
activation of slip events and the waiting times in between consecutive events as a superposition 
of Poisson processes. There are several possibilities to implement the KMC. The class @ref 
mepls::dynamics::KMC implements the method described in @cite DFCastellanos_CRP @cite Castellanos2019
 @cite Castellanos2018 @cite FernandezCastellanos2019.

We consider an extra ingredient in the model of thermal deformation. As we saw in previous 
tutorials, after a thermally activated slip event occurs, other slip systems might become 
unstable. We consider that mechanically unstable slip systems become active much faster than the 
time necessary for a new thermal activation somewhere in the whole system (this approximation 
breaks down if the temperature is too high, specifically if the material approaches a point such 
as the glass transition in glasses). Therefore, before triggering a new thermal event 
with the KMC, we apply an athermal relaxation as in the AQS simulations by calling @ref 
mepls::dynamics::relaxation. 



## Setting the temperature{#setting_temperature}

The model's temperature value is the rescaled one, \f$ T^{\prime} = k_{\rm B}T/V_{\rm a} \f$ where
 \f$ V_{\rm a} \f$ is of the same order as the product of the typical local plastic strain \f$ 
\Delta\gamma^{*}_{\rm pl} \f$ induced by an event and the volume occupied by the event. 
In this case, since we are working in 2D, it's an area. Let's denote by \f$ l^{*} \f$ the linear 
length of the event's region. The length of a mesoscale elements is \f$ l \f$, and at that 
scale the plastic strain is \f$ \Delta\gamma_{\rm pl} \f$. All quantities are related as

\f[  \Delta\gamma^{*}_{\rm pl} \cdot (l^{*})^2 = \Delta\gamma_{\rm pl} \cdot  l^2 \f]

If we use mesoscale elements that have the same size as the local plastic event, then \f$ l^{*}=l 
\f$, and therefore \f$ \Delta\gamma^{*}_{\rm pl} = \Delta\gamma_{\rm pl} \f$. The value of
\f$ \Delta\gamma_{\rm pl} \f$ is known since it is given by the input parameter `gamma`. 
Thus, if we set e.g. `gamma=0.05`, then \f$ \Delta\gamma^{*}_{\rm pl} = 0.05\f$. Now, let's say 
that the scale of a local plastic rearrangement is, for the microstructure that we are modelling, 1
 nm. In this case we have \f$ V_{\rm a} = \Delta\gamma^{*}_{\rm pl} \cdot (l^{*})^2 = 0.05 
 \textrm{nm}^2  \f$. The rescaled temperature is then

\f[  T^{\prime} = 3 \cdot 10^{-4} T \f]

which for a material undergoing thermal activation at, let's say \f$ T = 360 K \f$, is \f$ 
T^{\prime} = 0.1\f$.

As discussed above, this value denotes the typical amplitude of the local shear stress fluctuations 
induced by temperature. This value is a tenth of the typical slip threshold, given by
a scale parameter `lambda = 1.0`. Thus, thermal fluctuations are small compared to the 
thresholds, which means that the material is far away from a liquid-like behavior (i.e., its melting 
point). This fact justifies the applicability of the current modeling approach, based on solid 
mechanics, at that temperature. 
  
   

# The commented program{#comented_program}

We add two new headers. The header `omp.h` has functions 
that will provide us with information such as, e.g., the number running threads and the current 
thread id. The header `deal.II/base/conditional_ostream.h` defines an output stream, which will 
only accept the input from a specific thread.

```cpp
#include <example.h>
#include <mepls/utils.h>
#include <mepls/solver.h>
#include <mepls/system.h>
#include <mepls/history.h>
#include <mepls/dynamics.h>
#include <cmdparser.hpp>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/path_search.h>
#include <boost/property_tree/json_parser.hpp>

// to use openMP functions
#include <omp.h>

// to use a conditional output stream 
#include <deal.II/base/conditional_ostream.h>
```

The parameters struct is almost the same as in the previous tutorial step. We add the 
parameters of the creep test, namely the external stress `ext_stress_creep`, the duration of 
the test `time_limit_creep`, and the `temperature`. On the other hand, we have removed the parameter 
`filename`. The reason is that now the output filenames will be generated automatically by the 
program, based on the simulation parameters. Also, The parameter `verbose` allows us to switch on 
or off printing to the screen.  The parameter `n_rep` defines how many repetitions (each with 
a different random seed) of the simulation we want to run.

```cpp

struct Parameters
{
   unsigned int seed = 1234567;
   unsigned int n_rep = 1;
   unsigned int Nx = 32;
   unsigned int Ny = 32;
   double G = 30.;
   double nu = 0.3;
   double gamma = 0.05;
   double k = 6.;
   double strain_limit_aqs = 0.05;
   double time_limit_creep = 1e5;
   double ext_stress_creep = 0.;
   double lambda = 1.;
   double temperature = 1.;
   bool verbose = true;

   void declare_entries(dealii::ParameterHandler &prm)
   {
      prm.enter_subsection("Section1");

      prm.declare_entry("seed", mepls::utils::str::to_string(seed), dealii::Patterns::Integer(0), "");
      prm.declare_entry("n_rep", mepls::utils::str::to_string(n_rep), dealii::Patterns::Integer(0), "");
      prm.declare_entry("Nx", mepls::utils::str::to_string(Nx), dealii::Patterns::Integer(0), "");
      prm.declare_entry("Ny", mepls::utils::str::to_string(Ny), dealii::Patterns::Integer(0), "");
      prm.declare_entry("G", mepls::utils::str::to_string(G), dealii::Patterns::Double(0.0), "");
      prm.declare_entry("nu", mepls::utils::str::to_string(nu), dealii::Patterns::Double(0.0), "");
      prm.declare_entry("gamma", mepls::utils::str::to_string(gamma), dealii::Patterns::Double(0.0), "");
      prm.declare_entry("strain_limit_aqs", mepls::utils::str::to_string(strain_limit_aqs), dealii::Patterns::Double(0.0), "");
      prm.declare_entry("time_limit_creep", mepls::utils::str::to_string(time_limit_creep), dealii::Patterns::Double(0.0), "");
      prm.declare_entry("ext_stress_creep", mepls::utils::str::to_string(ext_stress_creep), dealii::Patterns::Double(0.0), "");
      prm.declare_entry("temperature", mepls::utils::str::to_string(temperature), dealii::Patterns::Double(0.0), "");
      prm.declare_entry("lambda", mepls::utils::str::to_string(lambda), dealii::Patterns::Double(0.0), "");
      prm.declare_entry("k", mepls::utils::str::to_string(k), dealii::Patterns::Double(0.0), "");
     prm.declare_entry("verbose", mepls::utils::str::to_string(verbose), dealii::Patterns::Bool(), "");

      prm.leave_subsection();
   }

   void load_entries(dealii::ParameterHandler &prm)
   {
      prm.enter_subsection("Section1");

      seed = prm.get_integer("seed");
      n_rep = prm.get_integer("n_rep");
      Nx = prm.get_integer("Nx");
      Ny = prm.get_integer("Ny");
      G = prm.get_double("G");
      nu = prm.get_double("nu");
      gamma = prm.get_double("gamma");
      k = prm.get_double("k");
      strain_limit_aqs = prm.get_double("strain_limit_aqs");
      time_limit_creep = prm.get_double("time_limit_creep");
      ext_stress_creep = prm.get_double("ext_stress_creep");
      temperature = prm.get_double("temperature");
      lambda = prm.get_double("lambda");
      verbose = prm.get_bool("verbose");

      prm.leave_subsection();
   }

   void load_file(const std::string &filename)
   {
      dealii::ParameterHandler prm;
      declare_entries(prm);
      prm.parse_input(filename);
      load_entries(prm);
   }

   void generate_file(const std::string &filename)
   {
      std::ofstream outfile(filename);
      dealii::ParameterHandler prm;
      declare_entries(prm);
      prm.print_parameters(outfile, dealii::ParameterHandler::OutputStyle::Text);
   }

};
```

This time, we will write two different histories, one for the creep and one for the AQS test.
We introduce a new level in the data hierarchy, adding `creep` and `AQS` between 
`Data` and the data sets. Thus, we have, e.g., `Data.creep.plastic_events` and `Data.AQS
.plastic_events`. The name of the output file is generated based on the simulation parameters. 
When generating the filename, the most important parameter is the `seed` parameter used for 
initializing the random number generator at the beginning of the `run` function. Adding the seed 
to the filename ensures that even if we perform many repetitions of the same simulation, files will
 never be overwritten unless they correspond to exactly the same simulation. Apart from that, the
 data is written in the same manner as in the previous tutorial step. 
 
```cpp
template<int dim>
void write_data(const mepls::history::History<dim> &creep_history,
                const mepls::history::History<dim> &aqs_history,
                const Parameters &p)
{
   boost::property_tree::ptree data_tree;

   data_tree.put("Name", "Step5");
   data_tree.put("Description", "System undergoing creep deformation, and then driven in "
                        "athermal quasistatic shear");

   data_tree.put("Parameters.dim", 2);
   data_tree.put("Parameters.seed", p.seed);
   data_tree.put("Parameters.Nx", p.Nx);
   data_tree.put("Parameters.Ny", p.Ny);
   data_tree.put("Parameters.G", p.G);
   data_tree.put("Parameters.nu", p.nu);
   data_tree.put("Parameters.gamma", p.gamma);
   data_tree.put("Parameters.lambda", p.lambda);
   data_tree.put("Parameters.temperature", p.temperature);
   data_tree.put("Parameters.k", p.k);
   data_tree.put("Parameters.ext_stress_creep", p.ext_stress_creep);
   data_tree.put("Parameters.time_limit_creep", p.time_limit_creep);
   data_tree.put("Parameters.strain_limit_aqs", p.strain_limit_aqs);

       // -------- creep history ----------
   {
      std::ostringstream plastic_events_csv;
      plastic_events_csv << "index,element,eigenstrain_00,eigenstrain_11,eigenstrain_01\n";
      for(auto &row : creep_history.plastic)
        plastic_events_csv << row.index << "," << row.element << "," << row.eigenstrain_00 << ","
                      << row.eigenstrain_11 << "," << row.eigenstrain_01 << "\n";

      data_tree.put("Data.creep.plastic_events", plastic_events_csv.str());


      std::ostringstream driving_events_csv;
      driving_events_csv << "index,dtime,dext_stress,dtotal_strain\n";
      for(auto &row : creep_history.driving)
        driving_events_csv << row.index << "," << row.dtime << ","
                      << row.dext_stress << "," << row.dtotal_strain << "\n";

      data_tree.put("Data.creep.driving_events", driving_events_csv.str());


      std::ostringstream macro_evolution_csv;
      macro_evolution_csv << "index,ext_stress,total_strain,time,av_vm_stress,av_vm_plastic_strain\n";
      for(auto &row : creep_history.macro_evolution)
        macro_evolution_csv << row.index << "," << row.ext_stress << ","
                      << row.total_strain << "," << row.time << ","
                      << row.av_vm_stress << "," << row.av_vm_plastic_strain << "\n";

      data_tree.put("Data.creep.macro_evolution", macro_evolution_csv.str());
   }

   // -------- aqs history ----------
   {
      std::ostringstream plastic_events_csv;
      plastic_events_csv << "index,element,eigenstrain_00,eigenstrain_11,eigenstrain_01\n";
      for(auto &row : aqs_history.plastic)
        plastic_events_csv << row.index << "," << row.element << "," << row.eigenstrain_00 << ","
                      << row.eigenstrain_11 << "," << row.eigenstrain_01 << "\n";

      data_tree.put("Data.AQS.plastic_events", plastic_events_csv.str());


      std::ostringstream driving_events_csv;
      driving_events_csv << "index,dtime,dext_stress,dtotal_strain\n";
      for(auto &row : aqs_history.driving)
        driving_events_csv << row.index << "," << row.dtime << ","
                      << row.dext_stress << "," << row.dtotal_strain << "\n";

      data_tree.put("Data.AQS.driving_events", driving_events_csv.str());


      std::ostringstream macro_evolution_csv;
      macro_evolution_csv << "index,ext_stress,total_strain,time,av_vm_stress,av_vm_plastic_strain\n";
      for(auto &row : aqs_history.macro_evolution)
        macro_evolution_csv << row.index << "," << row.ext_stress << ","
                      << row.total_strain << "," << row.time << ","
                      << row.av_vm_stress << "," << row.av_vm_plastic_strain << "\n";

      data_tree.put("Data.AQS.macro_evolution", macro_evolution_csv.str());
   }


   // create a descriptive filename
       std::ostringstream filename;

   filename << "Nx_" << p.Nx
            << "+gamma_" << std::fixed << std::setprecision(2) << p.gamma
            << "+lambda_" << p.lambda
            << "+k_" << p.k
           << "+G_" << p.G
           << "+T_" << p.temperature
           << "+ext_stress_" << p.ext_stress_creep
          << "+seed_" << p.seed
          << ".json";

   std::ofstream output_file( filename.str() );
   boost::property_tree::json_parser::write_json(output_file, data_tree);
   output_file.close();
}
```

The `run` function takes as an argument the parameters struct and a specialized output stream 
implemented by the `deal.II`, that will be described later. First, we set up the elements, the 
solver, and the system, as done in the previous tutorials. However, this time, the solver is 
configured to apply a stress-controlled load. Since we are 
using the solver @ref mepls::elasticity_solver::LeesEdwards<dim>, the stress-controlled load 
value refers to the xy-component of the externally applied stress. We start with a 
stress-controlled load since, before the AQS, we will perform the creep simulation. Later, when 
we perform the AQS test, we will change the driving conditions to strain-controlled.

```cpp
void run(const Parameters &p, dealii::ConditionalOStream & cout)
{
   constexpr unsigned dim = 2;
   std::mt19937 generator(p.seed);

   dealii::SymmetricTensor<4, dim> C = mepls::utils::tensor::make_isotropic_stiffness<dim>(p.G, p.nu);

   mepls::element::Vector<dim> elements;

   for(double n = 0; n < p.Nx * p.Ny; ++n)
   {
      example::element::Scalar<dim>::Config conf;
      conf.number = n;
      conf.gamma = p.gamma;
      conf.lambda = p.lambda;
      conf.k = p.k;
      conf.temperature = p.temperature;

      auto element = new example::element::Scalar<dim>(conf, generator);
      element->C(C);

      elements.push_back(element);
   }

   mepls::elasticity_solver::LeesEdwards<dim> solver(p.Nx, p.Ny, mepls::elasticity_solver::ControlMode::stress);
   for(auto &element : elements)
      solver.set_elastic_properties(element->number(), element->C());
   solver.setup_and_assembly();

   mepls::element::calculate_local_stress_coefficients_central(elements, solver);
   mepls::element::calculate_ext_stress_coefficients(elements, solver);

   mepls::system::Standard<dim> system(elements, solver, generator);
```

## Creep simulation{#creep}

We create a history object specific for the creep test, named `creep_history`. Then, for performing a 
creep simulation, we apply certain external stress and keep it fixed. To do this, we use the 
function @ref mepls::dynamics::fixed_load_increment, which applies to the system the 
desired load. Since the solver is currently in stress-controlled mode, the total strain will 
change during the deformation process, but the applied stress will not.

@note applying a specific load value is a straightforward operation since all we need to do is to call
@ref mepls::elasticity_solver::Solver<dim>::add_load_increment. However, the advantage of using 
MEPLS built-in calls such as @ref mepls::dynamics::fixed_load_increment is that the event 
associated with the load application will be recorded. Moreover, the elastic fields will 
be computed and the elements' stress tensors updated.
 
```cpp
   mepls::history::History<dim> creep_history("creep_history");
   system.set_history(creep_history);
   creep_history.add_macro(system);

   // apply an external (stress) load of amplitude p.ext_stress_creep
   mepls::dynamics::fixed_load_increment(p.ext_stress_creep, system);
   creep_history.add_macro(system);
```

We start the creep simulation by applying the Kinetic Monte Carlo method, which is 
implemented in the class @ref mepls::dynamics::KMC<dim>. The procedure is analogous to the 
AQS dynamics of previous tutorials, since we apply a perturbation followed 
by a relaxation. Here, the perturbation is a thermally activated event instead of an 
external load increment. The simulation loop stops when the time limit is reached.

```cpp
   mepls::dynamics::KMC<dim> kmc;

   mepls::utils::ContinueSimulation continue_creep;

   while( continue_creep() )
   {
      cout << system.macrostate["time"] << " " << system.macrostate["total_strain"] <<" "
           << system.macrostate["ext_stress"] << std::endl;

      kmc(system);
      creep_history.add_macro(system);

      mepls::dynamics::relaxation(system, continue_creep);
      creep_history.add_macro(system);

     continue_creep(system.macrostate["time"] < p.time_limit_creep, "creep time limit reached");
   }

    cout << continue_creep << std::endl;
```

## Transition step and AQS simulation{#transition}

When the creep simulation ends, we remove the external load. Immediately after the removal, some 
slip systems might become unstable. Therefore, we call @ref mepls::dynamics::relaxation to 
perform the associated slip events.
   
```cpp
   // remove the external load 
    mepls::dynamics::fixed_load_increment(-p.ext_stress_creep, system);
   creep_history.add_macro(system);

   // relax possible unstable slip system after the load change
   mepls::dynamics::relaxation(system, continue_creep);
   creep_history.add_macro(system);
```

The current configuration of slip thresholds and local stress tensors defines the initial
system's configuration during the AQS simulation. However, since we consider the AQS a process 
separated from sample preparation, we clean the current macroscale record. In this way, we start measuring the macroscale strain from zero. 

Currently, the solver operates in stress-controlled mode. However, we want to perform the AQS 
under strain-controlled conditions. We can change that in the existing solver bypassing 
@ref mepls::elasticity_solver::ControlMode::strain to 
@ref mepls::elasticity_solver::LeesEdwards<dim>::set_control_mode. Since the driving mode have 
changed, the local stress change per unit load increment is now different. Thus, we need to 
call @ref mepls::element::calculate_ext_stress_coefficients to inform the elements about this.
 
```cpp
    // in the AQS, we start measuring the macroscale strain from zero
   system.macrostate.clear();       

   // we switch the driving mode to strain-controlled during AQS
    solver.set_control_mode(mepls::elasticity_solver::ControlMode::strain);

    // since the type of driving conditions have changed, the local stress change induced by a unit
    // load increment is different. We need to re-compute the ext_stress_coefficients
    mepls::element::calculate_ext_stress_coefficients(elements, solver);
```

Now we perform the AQS simulation. The only difference with the previous tutorial 
step is that now we use the function @ref mepls::dynamics::extremal_dynamics_step. This function 
implements the quasistatic limit as accurately as possible within this model. Specifically, it 
performs an external load increment such that the barrier of one, and only one, slip system 
becomes zero. Therefore, after the load increment the system is infinitesimally close to 
undergoing plastic deformation. This critical corresponds to the minimal possible increment, in 
contrast with the previous tutorial step, where @ref mepls::dynamics::finite_extremal_dynamics_step performed a 
very small load but the probability of triggering several events at once was not zero. 

@note For using @ref mepls::dynamics::extremal_dynamics_step, having called @ref 
mepls::element::calculate_ext_stress_coefficients before is mandatory.

```cpp
       mepls::history::History<dim> aqs_history("aqs_history");
   system.set_history(aqs_history);
       aqs_history.add_macro(system);

   mepls::utils::ContinueSimulation continue_AQS_simulation;
   while(continue_AQS_simulation())
   {
      cout << system.macrostate["total_strain"] << " " << system.macrostate["ext_stress"] << std::endl;

      mepls::dynamics::extremal_dynamics_step(system);
      aqs_history.add_macro(system);

      mepls::dynamics::relaxation(system, continue_AQS_simulation);
      aqs_history.add_macro(system);

      continue_AQS_simulation(system.macrostate["total_strain"] < p.strain_limit_aqs, "AQS strain limit reached");
   }

      cout << continue_AQS_simulation << std::endl;

   for(auto &element : elements)
      delete element;

   // write the simulation data, using its own dedicated function
   write_data(creep_history, aqs_history, p);
   
} // run
```

The `main` function works in the same way as the previous tutorial step, but now we will run 
different repetitions in parallel using thread parallelism provided by openMP. Different 
repetitions denote different simulation runs with the same parameters but with a different random 
seed. For performing several repetitions, we call `run` within a loop. To perform repetitions in 
parallel, we enclose the loop within an openMP parallel region, as defined by `#pragma omp parallel` 
(see the openMP documentation). Within the parallel region, the process will run using several 
threads. Since we want a total number of repetitions defined by the `n_rep` parameter, each 
thread must run only that number divided by the number of threads.  

However, care must be taken to ensure that each parallel repetition differs from the repetitions
performed by the other threads. For that, we need to ensure that the random generator at 
the beginning of the `run` function is never initialized with the same seed. To this end, we 
consider the seed in the input parameters file as a master seed. We use the master seed to 
initialize a master random generator, which is used to generate different seeds passed to 
the `run` function. As a master generator, we use `std::rand`. This generator is of low quality, but 
that is not important for us since it only generates seeds that will initialize a much better 
`std::mt19937` generator during the simulation runs.

The master generator must not be used simultaneously by different threads, otherwise the seeds 
might be repeated. Therefore, we enclose it in a critical region, defined by the `#pragma 
critical`. 
Within this region, only one thread is allowed to enter at once. Since most of the time is spent 
inside the `run` function, the presence of the critical region has no impact on the performance.

The output stream that is passed to the `run` function ican be switched on or off. It is useful 
in a multi-threaded environment, since we can activate it for only one thread to avoid polluting 
the screen with simultaneous output.

```cpp
int main(int argc, char *argv[])
{
   // Read the command line arguments. We define the -f 'filename' to pass the
   // path to the parameters file
   cli::Parser parser(argc, argv);
   parser.set_optional<std::string>("f", "file", "./default.prm", "Name of the input configuration file");
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
      std::cout << "Configuration file " << parser.get<std::string>("f") << " created" << std::endl;
      return 1;
   }



	unsigned int n_rep = p.n_rep;
	if(n_rep < omp_get_max_threads())
		n_rep = omp_get_max_threads();

	// initialize the master engine with the master seed
	std::srand(p.seed);

	#pragma omp parallel
	{
		unsigned int n_threads = omp_get_max_threads();
		unsigned int id = omp_get_thread_num();
		unsigned int rep_per_thread = int( n_rep / n_threads );

		dealii::ConditionalOStream cout(std::cout, id==0 and p.verbose);

		for(unsigned int n = 0; n < rep_per_thread; ++n)
		{
			Parameters p_thread = p;
			#pragma critical
			{
				// generate a seed for a simulation run
				p_thread.seed = std::rand();
			};

			run(p_thread, cout);
		}

	}

   return 0;
}
```

   

# Results{#results_5}

We run the program as in @ref Step4, which generates a template parameters file,

```sh
$ ./run_sim
Configuration file ./default.prm created
```

then, we modify the default parameters file with the following values, 

```
# Listing of Parameters
# ---------------------
subsection Section1
  set G            = 30
  set Nx           = 32
  set Ny           = 32
  set n_rep        = 20
  set verbose      = true
  set gamma        = 0.05
  set k            = 6.00
  set temperature  = 0.1
  set lambda       = 1.00
  set nu           = 0.30
  set seed         = 123456
  set strain_limit_aqs = 0.05
  set time_limit_creep = 1e4
  set ext_stress_creep = 0.02
end
```

The value for the temperature is discussed in the introduction. During the creep simulation, we 
apply external stress of `ext_stress_creep = 0.02`. This value is very small in comparision 
with the thresholds (of the order of `lambda = 1.00`) and the stress fluctuations (of the order 
of `T = 0.1`), so it acts only as a small bias in the system's deformation leading to non-zero net deformation.

We set a number of repetitions `n_rep=20`. To perform those repetitions in parallel, we 
can set the number of threads we want to use by setting the environment variable 
`OMP_NUM_THREADS`. Thus, for using 4 threads we do

```bash
export OMP_NUM_THREADS=4
```

Also, it helps to set the variable

```bash
export OMP_PROC_BIND=TRUE
```

which pins the threads to specific cores avoiding continous switches and increases the performance 
about 20%. After setting these variables and modifying the parameters file, we run the program,

```sh
$ ./run_sim
0 0.000333333 0.02
0.43904 0.000308919 0.02
0.998616 0.000284505 0.02
1.76175 0.000308919 0.02
2.23352 0.000284505 0.02
2.29004 0.000260091 0.02
3.54874 0.000284505 0.02
3.67998 0.000260091 0.02
4.27786 0.000284505 0.02
4.64152 0.000308919 0.02
6.29096 0.000333333 0.02
6.66255 0.000357747 0.02
...
```

The result of the creep simulation can be animated with the Python script `animate.py`, located in 
the `MEPLS/python` directory. See @ref Step4 for a detailed explanation. We can run it as

```bash
python animate.py Nx_32+gamma_0.05+lambda_1.00+k_6.00+G_30.00+T_0.10+ext_stress_0.02+seed_1303402199.json /Data/creep/plastic_events /Data/creep/driving_events time vm_plastic_strain dstrain --labels '$t/\Delta t_0$' '$\varepsilon_{\rm vm}$ (%)' '$\varepsilon_{\rm vm}(\vec{r})$' --rescale 1 100 1
```

(note that the location and name of the output JSON file might differ). The call produces the 
following animation,

<center><img src="step5_plastic_strain_time_curve.gif" width="40%"></center>

The plot on the left shows the average von Mises plastic deformation vs. time curve. After a transient creep regime, we see how the system enters the stationary creep characterized by a 
constant strain rate. On the right, we can see the spatial activity. We cannot observe the 
formation of any clear patterns because the external stress is very low, as explained above. In 
this case, thermal activation leads to very noisy activity. 

The final state of the creep simulation was used as the initial state for the AQS one. To show 
the results of the AQS simulation, we do 

```bash
python animate.py Nx_32+gamma_0.05+lambda_1.00+k_6.00+G_30.00+T_0.10+ext_stress_0.02+seed_1303402199.json /Data/AQS/plastic_events /Data/AQS/driving_events total_strain ext_stress dstrain --labels '$\varepsilon_{\rm xy}$ (%)' '$\Sigma_{\rm xy}$ (MPa)' '$\varepsilon_{\rm vm}(\vec{r})$' --rescale 100 1000 1
```

which produces the following animation

<center><img src="step5_stress_strain_curve.gif" width="40%"></center>

We see that result is very similar to @ref Step4, where we considered that the slip thresholds 
where initially higher than after they are renewed. Here, the explanation is the same, but the 
initial thresholds configuration is the result of the creep dynamics. After the creep 
simulation, the system has undergone a certain plastic deformation at non-zero temperature and 
applied stress. In these conditions, the slip thresholds have experienced a survival bias, by 
which lower thresholds tend to be renewed by higher ones. Due to this bias, the statistics of the
 existing thresholds differ from the Weibull pdf from where they are renewed. Specifically, they 
 have a higher average. This state with higher thresholds is used as the initial state for the AQS 
simulation, leading to a stress overshoot and localization into a permanent shear band. The 
amplitude of the stress overshoot will depend on the temperature, the applied stress and 
the duration of the creep simulation (but it if reaches the stationary state, it won't 
change anymore). Thus, although the results are similar to @ref Step4, here we have built a model 
much more powerful, which provides us with a great insight into the effects of a 
material's thermal and mechanical history on its brittle or ductile response.



# The complete program{#full_5}

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
#include <mepls/dynamics.h>
#include <cmdparser.hpp>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/path_search.h>
#include <boost/property_tree/json_parser.hpp>

// new headers
#include <omp.h>
#include <deal.II/base/conditional_ostream.h>


struct Parameters
{
   unsigned int seed = 1234567;
   unsigned int n_rep = 1;
   unsigned int Nx = 32;
   unsigned int Ny = 32;
   double G = 30.;
   double nu = 0.3;
   double gamma = 0.05;
   double k = 6.;
   double strain_limit_aqs = 0.05;
   double time_limit_creep = 1e5;
   double ext_stress_creep = 0.;
   double lambda = 1.;
   double temperature = 1.;
   bool verbose = true;

   void declare_entries(dealii::ParameterHandler &prm)
   {
      prm.enter_subsection("Section1");

      prm.declare_entry("seed", mepls::utils::str::to_string(seed), dealii::Patterns::Integer(0), "");
      prm.declare_entry("n_rep", mepls::utils::str::to_string(n_rep), dealii::Patterns::Integer(0), "");
      prm.declare_entry("Nx", mepls::utils::str::to_string(Nx), dealii::Patterns::Integer(0), "");
      prm.declare_entry("Ny", mepls::utils::str::to_string(Ny), dealii::Patterns::Integer(0), "");
      prm.declare_entry("G", mepls::utils::str::to_string(G), dealii::Patterns::Double(0.0), "");
      prm.declare_entry("nu", mepls::utils::str::to_string(nu), dealii::Patterns::Double(0.0), "");
      prm.declare_entry("gamma", mepls::utils::str::to_string(gamma), dealii::Patterns::Double(0.0), "");
      prm.declare_entry("strain_limit_aqs", mepls::utils::str::to_string(strain_limit_aqs), dealii::Patterns::Double(0.0), "");
      prm.declare_entry("time_limit_creep", mepls::utils::str::to_string(time_limit_creep), dealii::Patterns::Double(0.0), "");
      prm.declare_entry("ext_stress_creep", mepls::utils::str::to_string(ext_stress_creep), dealii::Patterns::Double(0.0), "");
      prm.declare_entry("temperature", mepls::utils::str::to_string(temperature), dealii::Patterns::Double(0.0), "");
      prm.declare_entry("lambda", mepls::utils::str::to_string(lambda), dealii::Patterns::Double(0.0), "");
      prm.declare_entry("k", mepls::utils::str::to_string(k), dealii::Patterns::Double(0.0), "");
	  prm.declare_entry("verbose", mepls::utils::str::to_string(verbose), dealii::Patterns::Bool(), "");

      prm.leave_subsection();
   }

   void load_entries(dealii::ParameterHandler &prm)
   {
      prm.enter_subsection("Section1");

      seed = prm.get_integer("seed");
      n_rep = prm.get_integer("n_rep");
      Nx = prm.get_integer("Nx");
      Ny = prm.get_integer("Ny");
      G = prm.get_double("G");
      nu = prm.get_double("nu");
      gamma = prm.get_double("gamma");
      k = prm.get_double("k");
      strain_limit_aqs = prm.get_double("strain_limit_aqs");
      time_limit_creep = prm.get_double("time_limit_creep");
      ext_stress_creep = prm.get_double("ext_stress_creep");
      temperature = prm.get_double("temperature");
      lambda = prm.get_double("lambda");
      verbose = prm.get_bool("verbose");

      prm.leave_subsection();
   }

   void load_file(const std::string &filename)
   {
      dealii::ParameterHandler prm;
      declare_entries(prm);
      prm.parse_input(filename);
      load_entries(prm);
   }

   void generate_file(const std::string &filename)
   {
      std::ofstream outfile(filename);
      dealii::ParameterHandler prm;
      declare_entries(prm);
      prm.print_parameters(outfile, dealii::ParameterHandler::OutputStyle::Text);
   }

};


template<int dim>
void write_data(const mepls::history::History<dim> &creep_history,
                const mepls::history::History<dim> &aqs_history,
                const Parameters &p)
{
   boost::property_tree::ptree data_tree;

   data_tree.put("Name", "Step5");
   data_tree.put("Description", "System undergoing creep deformation, and then driven in "
								"athermal quasistatic shear");

   data_tree.put("Parameters.dim", 2);
   data_tree.put("Parameters.seed", p.seed);
   data_tree.put("Parameters.Nx", p.Nx);
   data_tree.put("Parameters.Ny", p.Ny);
   data_tree.put("Parameters.G", p.G);
   data_tree.put("Parameters.nu", p.nu);
   data_tree.put("Parameters.gamma", p.gamma);
   data_tree.put("Parameters.lambda", p.lambda);
   data_tree.put("Parameters.temperature", p.temperature);
   data_tree.put("Parameters.k", p.k);
   data_tree.put("Parameters.ext_stress_creep", p.ext_stress_creep);
   data_tree.put("Parameters.time_limit_creep", p.time_limit_creep);
   data_tree.put("Parameters.strain_limit_aqs", p.strain_limit_aqs);

   	// -------- creep history ----------
	{
	   std::ostringstream plastic_events_csv;
	   plastic_events_csv << "index,element,eigenstrain_00,eigenstrain_11,eigenstrain_01\n";
	   for(auto &row : creep_history.plastic)
		  plastic_events_csv << row.index << "," << row.element << "," << row.eigenstrain_00 << ","
							 << row.eigenstrain_11 << "," << row.eigenstrain_01 << "\n";

	   data_tree.put("Data.creep.plastic_events", plastic_events_csv.str());


	   std::ostringstream driving_events_csv;
	   driving_events_csv << "index,dtime,dext_stress,dtotal_strain\n";
	   for(auto &row : creep_history.driving)
		  driving_events_csv << row.index << "," << row.dtime << ","
							 << row.dext_stress << "," << row.dtotal_strain << "\n";

	   data_tree.put("Data.creep.driving_events", driving_events_csv.str());


	   std::ostringstream macro_evolution_csv;
	   macro_evolution_csv << "index,ext_stress,total_strain,time,av_vm_stress,av_vm_plastic_strain\n";
	   for(auto &row : creep_history.macro_evolution)
		  macro_evolution_csv << row.index << "," << row.ext_stress << ","
							 << row.total_strain << "," << row.time << ","
							 << row.av_vm_stress << "," << row.av_vm_plastic_strain << "\n";

	   data_tree.put("Data.creep.macro_evolution", macro_evolution_csv.str());
	}

	// -------- aqs history ----------
	{
	   std::ostringstream plastic_events_csv;
	   plastic_events_csv << "index,element,eigenstrain_00,eigenstrain_11,eigenstrain_01\n";
	   for(auto &row : aqs_history.plastic)
		  plastic_events_csv << row.index << "," << row.element << "," << row.eigenstrain_00 << ","
							 << row.eigenstrain_11 << "," << row.eigenstrain_01 << "\n";

	   data_tree.put("Data.AQS.plastic_events", plastic_events_csv.str());


	   std::ostringstream driving_events_csv;
	   driving_events_csv << "index,dtime,dext_stress,dtotal_strain\n";
	   for(auto &row : aqs_history.driving)
		  driving_events_csv << row.index << "," << row.dtime << ","
							 << row.dext_stress << "," << row.dtotal_strain << "\n";

	   data_tree.put("Data.AQS.driving_events", driving_events_csv.str());


	   std::ostringstream macro_evolution_csv;
	   macro_evolution_csv << "index,ext_stress,total_strain,time,av_vm_stress,av_vm_plastic_strain\n";
	   for(auto &row : aqs_history.macro_evolution)
		  macro_evolution_csv << row.index << "," << row.ext_stress << ","
							 << row.total_strain << "," << row.time << ","
							 << row.av_vm_stress << "," << row.av_vm_plastic_strain << "\n";

	   data_tree.put("Data.AQS.macro_evolution", macro_evolution_csv.str());
	}


	// create a descriptive filename
   	std::ostringstream filename;

	filename << "Nx_" << p.Nx
	         << "+gamma_" << std::fixed << std::setprecision(2) << p.gamma
	         << "+lambda_" << p.lambda
	         << "+k_" << p.k
		     << "+G_" << p.G
		     << "+T_" << p.temperature
		     << "+ext_stress_" << p.ext_stress_creep
			 << "+seed_" << p.seed
			 << ".json";

   std::ofstream output_file( filename.str() );
   boost::property_tree::json_parser::write_json(output_file, data_tree);
   output_file.close();
}


void run(const Parameters &p, dealii::ConditionalOStream & cout)
{
   //----- SET UP ------

   constexpr unsigned dim = 2;
   std::mt19937 generator(p.seed);

   dealii::SymmetricTensor<4, dim> C = mepls::utils::tensor::make_isotropic_stiffness<dim>(p.G, p.nu);

   mepls::element::Vector<dim> elements;

   for(double n = 0; n < p.Nx * p.Ny; ++n)
   {
      example::element::Scalar<dim>::Config conf;
      conf.number = n;
      conf.gamma = p.gamma;
      conf.lambda = p.lambda;
      conf.k = p.k;
      conf.temperature = p.temperature;

      auto element = new example::element::Scalar<dim>(conf, generator);
      element->C(C);

      elements.push_back(element);
   }

   mepls::elasticity_solver::LeesEdwards<dim> solver(p.Nx, p.Ny, mepls::elasticity_solver::ControlMode::stress);
   for(auto &element : elements)
      solver.set_elastic_properties(element->number(), element->C());
   solver.setup_and_assembly();

   mepls::element::calculate_local_stress_coefficients_central(elements, solver);
   mepls::element::calculate_ext_stress_coefficients(elements, solver);

   mepls::system::Standard<dim> system(elements, solver, generator);


   //----- CREEP DEFORMATION ------

   mepls::history::History<dim> creep_history("creep_history");
   system.set_history(creep_history);
   creep_history.add_macro(system);

   // apply an external (stress) load of amplitude p.ext_stress_creep
   mepls::dynamics::fixed_load_increment(p.ext_stress_creep, system);
   creep_history.add_macro(system);

   mepls::dynamics::KMC<dim> kmc;

   mepls::utils::ContinueSimulation continue_creep;

	while( continue_creep() )
	{
		cout << system.macrostate["time"] << " " << system.macrostate["total_strain"] <<" "
			  << system.macrostate["ext_stress"] << std::endl;

		kmc(system);
		creep_history.add_macro(system);

		mepls::dynamics::relaxation(system, continue_creep);
		creep_history.add_macro(system);

	  continue_creep(system.macrostate["time"] < p.time_limit_creep, "creep time limit reached");
	}

	 cout << continue_creep << std::endl;


   //----- TRANSITION TO AQS ------

	// remove the external load
    mepls::dynamics::fixed_load_increment(-p.ext_stress_creep, system);
	creep_history.add_macro(system);

	// relax possible unstable slip system after the load change
	mepls::dynamics::relaxation(system, continue_creep);
	creep_history.add_macro(system);

    // in the AQS, we start measuring the macroscale strain from zero
	system.macrostate.clear();

	for(auto & element : elements)
		element->state_to_prestress();

	// we switch the driving mode to strain-controlled during AQS
    solver.set_control_mode(mepls::elasticity_solver::ControlMode::strain);

    // since the type of driving conditions have changed, the local stress change induced by a unit
    // load increment is different. We need to re-compute the ext_stress_coefficients
    mepls::element::calculate_ext_stress_coefficients(elements, solver);


    //----- ATHERMAL QUASISTATIC SHEAR ------

   	mepls::history::History<dim> aqs_history("aqs_history");
  	system.set_history(aqs_history);
   	aqs_history.add_macro(system);

   mepls::utils::ContinueSimulation continue_AQS_simulation;
   while(continue_AQS_simulation())
   {
		cout << system.macrostate["total_strain"] << " " << system.macrostate["ext_stress"] << std::endl;

      mepls::dynamics::extremal_dynamics_step(system);
      aqs_history.add_macro(system);

      mepls::dynamics::relaxation(system, continue_AQS_simulation);
      aqs_history.add_macro(system);

      continue_AQS_simulation(system.macrostate["total_strain"] < p.strain_limit_aqs, "AQS strain limit reached");
   }

	   cout << continue_AQS_simulation << std::endl;

   for(auto &element : elements)
      delete element;

	// write the simulation data, using its own dedicated function
	write_data(creep_history, aqs_history, p);
}


int main(int argc, char *argv[])
{
   // Read the command line arguments. We define the -f 'filename' to pass the
   // path to the parameters file
   cli::Parser parser(argc, argv);
   parser.set_optional<std::string>("f", "file", "./default.prm", "Name of the input configuration file");
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
      std::cout << "Configuration file " << parser.get<std::string>("f") << " created" << std::endl;
      return 1;
   }



	unsigned int n_rep = p.n_rep;
	if(n_rep < omp_get_max_threads())
		n_rep = omp_get_max_threads();

	// initialize the master engine with the master seed
	std::srand(p.seed);

	#pragma omp parallel
	{
		unsigned int n_threads = omp_get_max_threads();
		unsigned int id = omp_get_thread_num();
		unsigned int rep_per_thread = int( n_rep / n_threads );

		dealii::ConditionalOStream cout(std::cout, id==0 and p.verbose);

		for(unsigned int n = 0; n < rep_per_thread; ++n)
		{
			Parameters p_thread = p;
			#pragma critical
			{
				// generate a seed for a simulation run
				p_thread.seed = std::rand();
			};

			run(p_thread, cout);
		}

	}

   return 0;
}
```
