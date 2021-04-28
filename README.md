# MEPLS
[https://mepls.davidfcastellanos.com/](https://mepls.davidfcastellanos.com/)


## Table of contents

* [About MEPLS](#about-MEPLS)
* [Features](#features)
* [Project structure](#project-structure)
* [How to build](#how-to-build)
* [Documentation](#documentation)
* [Contact](#contact)
* [License](#license)


## About MEPLS
MEPLS (Mesoscale Elasto-Plasticity Simulator) is an open-source C++ framework 
for simulating the plastic deformation of materials based on the combination of 
stochastic processes and solid mechanics. Its goal is to provide modular, 
efficient, and easy-to-use tools to build simulations for a broad range of
physical scenarios. MEPLS abstracts away unnecessary complexity, allowing the 
user to focus on the physics of the model. Whenever possible, it uses existing 
well-established software, such as the [deal.II] library for the Finite Element 
Method. To achieve these goals, MEPLS is developed around three pillars:

  * **Efficient** --- it is fully implemented in C++, with a back-end based on 
    well-established efficient libraries such as deal.II and UMFPACK.

  * **User-friendly** --- it is easy to build and has minimal dependencies. It 
    is well-documented and has a gallery of examples.

  * **Extensible** --- it is object-oriented and template-based. Its modularity 
    and simplicity facilitate the task of extending it according to the user's 
    needs.


##### Background & motivation
At a range of scales between the micro and the macroscale, a material's plastic 
deformation can be described with great independence of microscopic details at 
the so-called mesocale. Models operating at the mesoscale are an excellent tool 
for establishing links between micro and macroscale plastic deformation.
Mesoscale models divide a material into mesoscale subdomains, known as elements.
The state of the elements is defined by internal variables in terms of continuum
mechanics. Their behavior is defined according to a set of simple local rules,
representing only the most relevant characteristics of the plastic deformation 
mechanisms and neglecting fine details. Such a simplified description means a 
loss of information and imperfect knowledge about the elements' state. To 
represent this uncertainty, the rules governing the behavior of the elements are
stochastic processes. Finally, elements influence each other's behavior due to 
elastic fields. This interaction leads to a highly correlated system that
exhibits complex non-linear behavior. In these situations, analytical solutions 
are not available, so we must resort to numerical methods.
 
MEPLS is a framework that provides easy-to-use numerical tools to build mesoscale
models of the plastic deformation of materials under a broad range of physical 
scenarios. An alternative description is that MEPLS allows the statistical 
sampling of a mesoscale model's configuration space with Markov chain Monte Carlo
methods. MEPLS enables the user to study the impact that microstructural 
properties have on the macroscale elasto-plastic behavior of a material. This 
question is especially challenging to answer when those properties are spatially
heterogeneous and statistically distributed, which is, in turn, the most common
situation. The built-in MEPLS tools abstract away unnecessary complexity 
associated with model specification, allowing the user to focus on the physics 
of the problem. Thus, the user can easily study the influence of: the local 
statistics of elastic constants, yield criteria, plastic deformation, structural
evolution, etc.; the driving protocol: strain-controlled tests, creep, cyclic 
loading, etc.; the loading mode: simple shear, compression, bending, etc. 
Moreover, suppose a feature is not readily available. In that case, chances are
that the user can seamlessly implement it thanks to MEPLS extensibility, 
availability of examples, and thorough documentation.

You can find more information about MEPLS on the [project's website]. 


## Features
MEPLS has some built-in features, while others have been implemented in some of 
the example models. If you are looking for some of these features, you are lucky
since you can use them out-of-the-box. Otherwise, they might be implemented in 
future releases, or you can do it yourself. Some of them might fall into the 
category for which MEPLS extensibility was planned and are therefore very 
straightforward to implement. Others may imply a significant amount of effort. 
Here is a list of currently available and unavailable features. If a feature is 
unavailable, you can see an estimation of how difficult it would be to implement
it. Of course, since the list of imaginable features is almost infinite, the 
currently unavailable features will always outnumber the available ones.  

#### Available
  * Control of the eigenstrain field

  * Control over microstructural properties and their evolution

  * Fully tensorial stress and strain

  * 2D simulations with quadrilateral structured meshes

  * Periodic boundary conditions

  * Quasistatic strain-driven tests

  * Kinetic Monte Carlo, Metropolis-Hastings, and extremal dynamics for the 
    activation of plastic events

  * Relaxation tests by thermally-activated plastic activity

  * Quenched heterogeneous elastic properties

  * Local plastic yield criteria: von Mises, Drucker-Prager, and slip-based 
    crystal plasticity

  * Mean-field elastic & randomized elastic interaction kernels

#### Available soon
  * Different external loading conditions

  * Quasistatic stress-driven tests

  * Creep by thermally-activated plastic activity

#### Currently unavailable but: 

##### Easy to implement
  * Any other local plastic yield criterion such as Mohr-Coulomb

  * Any law for the evolution of local plastic properties: strain-softening, 
    local densification, etc.

  * Fracture by strain-softening or by a reduction of local stiffness

  * Cyclic loading

  * Different protocols to drive the material, or combinations of them (e.g., 
    several cycles of loading followed by a stress relaxation, or creep followed
    by a strain-driven test)

  * Creating material sub-domains with different properties to simulate, e.g.,
    multi-phase composites

  * Evolution of homogeneous elastic properties with time or strain

##### Challenging to implement
  * Non-quasistatic strain- or stress-driven tests

  * Evolution of non-homogeneous elastic properties (computationally expensive 
    and would probably require changing the FEM direct solver for an iterative 
    one provided by the deal.II back-end)

  * Viscosity effects

  * 3D 

##### Difficult to implement
  * Non-structured meshes

  * Triangular meshes (impossible with the deal.II library, but adding a 
    different FEM library could make it possible; MEPLS would still depend on 
    deal.II tensor classes to work)

  * Internal parallelization of individual simulation runs


##### Implementation affecting only the FEM solver
The following features are easily compatible with MEPLS, and their implementation
deals mainly with creating a new deal.II-based FEM solver specialized for the 
task. Therefore, the difficulty of implementing them depends on the user's 
knowledge of the deal.II library (note: these features would also make simulation
runs very computationally demanding, difficulting the statistical sampling of 
the models):

  * Non-linear elasticity

  * Surface contact problem

  * Elastic waves propagation

  * Finite strains

  * Material advection


## Project structure
MEPLS is the combination of a header-only library and a set of bundled
example models illustrating how to apply and combine different classes and 
algorithms provided by the headers. If none of the existing models is suitable 
for a user's needs, a new model can be easily implemented following the example 
models and the project's documentation. This task will involve implementing some
classes derived from the abstract classes provided by MEPLS. The implementation
of the derived classes focuses on the physics of the model that the user is 
interested in. Once this is done, the user can rely on MEPLS to handle the rest.

#### MEPLS core
The core of MEPLS is the header-only library located in 
`/path/to/MEPLS_dir/inc/MEPLS`. The available tools are defined in the namespace
`mepls`. Within that namespace, there are several namespaces containing classes 
and functions with well-differentiated reponsabilities:

| Namespace | Responsibility of the members |
| ------ | ------ |
| `element` | represent the material's mesoscale subdomains, a.k.a elements |
| `slip` | implement the physics of slip activation within the elements |
| `elasticity_solver` | compute the elastic fields, i.e., how the elements influence each other |
| `history` | record, continuously in time, the evolution and changes in the system |
| `system` | manage the elements, the solver, and the history  |
| `dynamics` | algorithms for the dynamical evolution of the system |
| `snapshot` | take snapshots of the elements' state at certain discrete time steps |
| `patches` | sample the local elasto-plastic response of the system at different length scales |

#### Existing models & key features
The bundled models are located in `/path/to/MEPLS_dir/models/`. The available 
models are:

  * ESPCI:

    * Random local slip planes as glass microstructure
    
    * Monte Carlo simulation of a super-cooled liquid
    
    * Athermal quasistatic strain-controlled test
    
    * Unloading and reloading
    
    * Sampling local elasto-plastic response using the patch method


## How to build

####Linux
MEPLS uses CMake for building. The only MEPLS dependency is the [deal.II] library,
with version between 8.5 and 9.1. The recommended version is the 9.0, which can 
be found in this repository. A minimal build of deal.II is enough. To build it, 
unpack the deal.II sources in some directory `/path/to/dealii/sources`. Then, to 
configure and build it in, e.g., `/path/to/dealii/build`, you can do:

```sh
mkdir /path/to/dealii/build
cd /path/to/dealii/build
cmake /path/to/dealii/sources
make # use the argument -j<N> for a build using N parallel threads
make test
```

Now, we can compile some of the MEPLS models. First, unpack MEPLS in some 
directory `/path/to/MEPLS`. Let’s build the model located in 
`/path/to/MEPLS/models/example`. To build it in `/path/to/model/build`, simply
do:

```sh
mkdir /path/to/model/build
cd /path/to/model/build
cmake /path/to/MEPLS/models/example -DDEAL_II_DIR=/path/to/dealii/build
make
```

After that, an executable named `run_sim` should appear in `/path/to/model/build`.

Each model is built using its own `CMakeLists.txt`. Thus, in general, cmake 
configuration flags are model-specific and are defined by the model’s author. 
However, there are two flags that control MEPLS built-in capabilities, namely:

  * `DDEFINE_DEBUG=[on/off]` defines whether to build the model using the debug
    mode. In debug mode, MEPLS performs intensive checks, which can reduce the 
    performance significantly.
    
  * `DOPENMP=[on/off]` will allow performing several simulation runs in parallel
    using OpenMP.

Details on how to control the parameters of the parallelization can be found in 
the [OpenMP](https://www.openmp.org/) documentation. The most important parameter
is the number of threads, which is specified by the environment variable 
`OMP_NUM_THREADS`. It is recommended that we also set the variable 
`OMP_PROC_BIND=TRUE` to pin the threads to specific CPU cores and prevent 
continous switches, which can introduce considerable overhead. Currently, MEPLS
does not support internal parallelization. Therefore, a single simulation run 
cannot run in parallel.

####Windows
Building on Windows has not been tested yet. You can follow [deal.II]'s own building
instructions for Windows, after which building a MEPLS model should be trivial.

####macOS
Building on macOS has not been tested yet. You can follow [deal.II]'s own building
instructions for macOS, after which building a MEPLS model should be trivial.


## Documentation
You can access the documentation online on the [project's website]. 
Alternatively, you can generate the documentation in different formats using 
[Doxygen](https://www.doxygen.nl/index.html) and the generation script located 
in `/path/to/MEPLS_dir/inc/MEPLS/doc`.


## Contact
You can contact the author for, e.g., reporting issues and bugs or submitting 
requests and ideas through the contact methods provided on the 
[author's website].


## License
MEPLS is open source. You can freely use it, redistribute it, and/or modify it
under the terms of the Creative Commons Attribution 4.0 International Public 
License. The full text of the license can be found in the file LICENSE at the 
top level of the MEPLS distribution.
 
Copyright (C) 2020  - David Fernández Castellanos.


   [deal.II]: <https://www.dealii.org/>
   [project's website]: <https://mepls.davidfcastellanos.com/>
   [author's website]: <https://www.davidfcastellanos.com/>
   
____