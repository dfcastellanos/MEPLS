# MEPLS

## Table of contents

* [About](#about)
* [Features](#features)
* [Documentation and tutorial](#documentation-and-tutorial)
* [How to build](#how-to-build)
* [Contact](#contact)
* [License](#license)


## About
MEPLS (Mesoscale Elasto-Plasticity Simulator) is an open-source C++ framework 
for simulating the plastic deformation of materials based on the combination of 
stochastic processes and solid mechanics. Its goal is to provide modular, 
efficient, and easy-to-use tools to build simulations for a broad range of
physical scenarios. MEPLS abstracts away unnecessary complexity, allowing the 
user to focus on the physics of the model. Whenever possible, it uses existing 
well-established software, such as the [deal.II] library for the Finite Element 
Method. To achieve these goals, MEPLS is developed around three pillars:

  * **Efficient** — it is fully implemented in C++, with a back-end based on 
    well-established efficient libraries such as deal.II and UMFPACK.

  * **User-friendly** — it is easy to build and has minimal dependencies. It 
    is well-documented and has a gallery of examples.

  * **Extensible** — it is object-oriented and template-based. Its modularity 
    and simplicity facilitate the task of extending it according to the user's 
    needs.

MEPLS provides easy-to-use tools to build simulations that statistically sample mesoscale models of 
the plastic deformation of materials. The built-in MEPLS tools abstract away unnecessary 
complexity associated with the model implementation, allowing the user to focus on the physics of
the problem.

MEPLS enables the user to study the impact that microstructural properties have on the macroscale
elasto-plastic behavior of a material. This question is especially challenging to answer when 
those properties are spatially heterogeneous and statistically distributed, which is, in turn, 
the most common situation. MEPLS allows to study the influence that many different model features
 or ingredients may have on a material's mechanical properties. Example of such features are:

* the statistics of local elastic properties, local yield criteria, plastic deformation, structural 
evolution, etc.
* the driving protocol: strain-controlled tests, creep, cyclic loading, etc.
* the loading mode: simple shear, compression, bending, etc.

You can find an extensive list in the features section.


MEPLS is developed and maintained by David Fernández Castellanos. You can find more information 
about MEPLS on the [project's website]. 

## Features
MEPLS has some built-in features, while others have been implemented in some of 
the bundled tutorial and gallery models. If you are looking for some of these features, you are 
lucky since you can use them out-of-the-box. Otherwise, they might be implemented in 
future releases, or you can do it yourself (to do it yourself, you can get familiar with MEPLS
by following its tutorial, taking a look at the gallery of models, and reading its official 
documentation). Some of the desired features might fall into the category for which MEPLS 
extensibility was planned and are straightforward to implement. Others may imply a
significant amount of effort. Here is a list of currently available and unavailable features. If 
a feature is unavailable, you can see an estimation of how difficult it would be to implement
it.


### Available
  * Manipulation of the plastic strain field and the microstructural properties

  * Fully tensorial elastic fields

  * 2D simulations with quadrilateral structured meshes

  * Periodic boundary conditions
  
  * Different algorithms for the activation of plastic slips such as, e.g.,
    Kinetic Monte Carlo, Metropolis-Hastings or extremal dynamics

  * Quasistatic strain- and stress-controlled tests

  * Creep by thermally-activated plastic activity
  
  * Relaxation tests by thermally-activated plastic activity

  * Quenched heterogeneous elastic properties

  * Local plastic yield criteria such as, e.g., von Mises, Drucker-Prager or slip-based 
    crystal plasticity

  * Mean-field elastic & randomized elastic interaction kernels

### Available soon
  * Different external loading conditions


### Currently unavailable: 

#### Easy to implement
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

##### The implementation deals only with the FEM solver
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


## Documentation and tutorial

You can access the documentation and the tutorial online on the [project's website]. 
Alternatively, you can generate them using [Doxygen](https://www.doxygen.nl/index.html) and the 
generation script located in `/path/to/MEPLS_dir/doc`.


## How to build
 
MEPLS is a header-only library. However, the distribution includes a set of bundled models that compose the
tutorial and the gallery. All those models can be built according to the following instructions. 
The instructions can also serve as a guide to build your own model.

First of all, MEPLS depends on the [deal.II] library for using the Finite Element Method. 
Compatible deal.II versions are between the 8.5 and the 9.1. The recommended is the 9.0, which 
can be found in [this repository](https://github.com/dealii/dealii/tree/dealii-9.0). You will need
the LAPACK library, which is very standard, and you might already have it installed in your 
system. If not, you can install it from your OS distribution repositories (for example, for 
Ubuntu-based systems, do `sudo apt-get install liblapack-dev`). 

Deal.II is built using [CMake](https://cmake.org/). MEPLS bundled models will also use it.  

#### Linux  

To build deal.II, unpack its sources in some directory `/path/to/dealii/sources`. Then, to 
configure and build it in, e.g., `/path/to/dealii/build`, you can do:

```sh
mkdir /path/to/dealii/build
cd /path/to/dealii/build
cmake /path/to/dealii/sources -DDEAL_II_WITH_UMFPACK=on
make # use the argument -j<N> for a build using N parallel threads
```
Alternatively, you can do `make install` to install the library in you system, but is not 
mandatory.

If you have any doubts about the process of building deal.II, you can find a more extensive 
explanation on its documentation, see [deal.II installation](https://www.dealii.org/current/readme.html#installation).

Now, you can compile some of the MEPLS models. First, unpack MEPLS in some 
directory `/path/to/MEPLS_dir`. The tutorial steps are located int `/path/to/MEPLS_dir/tutorial` 
and the gallery models in `/path/to/MEPLS_dir/gallery`. Let’s build the model located in 
`/path/to/MEPLS_dir/tutorial/step5`. To build it in `/path/to/model/build`, simply
do:

```sh
mkdir /path/to/model/build
cd /path/to/model/build
cmake /path/to/MEPLS_dir/tutorial/step5 -DDEAL_II_DIR=/path/to/dealii/build -DMEPLS_DIR=/path/to/MEPLS_dir
make
```

Instead of adding the flags `-DDEAL_II_DIR` and `-DMEPLS_DIR`, you can also set the 
corresponding environment variables as `export DEAL_II_DIR=/path/to/dealii/build` and `export 
MEPLS_DIR=/path/to/MEPLS_dir`. Also, if you chose to install deal.II in your system, CMake should
find it automatically.

After calling make, an executable named `run_sim` should appear in `/path/to/model/build`.

Each model is built using its own `CMakeLists.txt`, so different models might accept different 
configuration flags. However, the flag `DEFINE_DEBUG` turns on MEPLS debug mode,

  * `-DDEFINE_DEBUG=[on/off]` defines whether to build the model using the debug
    mode. In debug mode, MEPLS performs intensive checks, which can reduce the 
    performance significantly. If the flag is not set, the default value is `off`.
    

#### Windows
Building on Windows has not been tested yet. You can follow [deal.II]'s building
instructions for Windows, after which building a MEPLS model should be straightforward.

#### macOS
Building on macOS has not been tested yet. You can follow [deal.II]'s building
instructions for macOS, after which building a MEPLS model should be straightforward.


## Contact
MEPLS is developed and maintained by David Fernández Castellanos. You can report issues and bugs 
in the [project's repository](https://github.com/dfcastellanos/MEPLS). You can contact the author 
through the methods provided on the [author's website] for longer discussions regarding, e.g., 
requests and ideas for the project or if you need some help to use MEPLS.


## License
MEPLS is open source. You can freely use it, redistribute it, and/or modify it
under the terms of the Creative Commons Attribution 4.0 International Public 
License. The full text of the license can be found in the file LICENSE at the 
top level of the MEPLS distribution.
 
Copyright (C) 2020  - David Fernández Castellanos.


   [deal.II]: <https://www.dealii.org/>
   [project's website]: <https://mepls.dfcastellanos.com/>
   [author's website]: <https://www.dfcastellanos.com/contact>
   
____
