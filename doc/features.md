@page Features

# Features
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



## Available
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

## Available soon
  * Different external loading conditions

  * Quasistatic stress-driven tests

  * Creep by thermally-activated plastic activity



## Currently unavailable but: 

### Easy to implement
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

### Challenging to implement
  * Non-quasistatic strain- or stress-driven tests

  * Evolution of non-homogeneous elastic properties (computationally expensive 
    and would probably require changing the FEM direct solver for an iterative 
    one provided by the deal.II back-end)

  * Viscosity effects

  * 3D 

### Difficult to implement
  * Non-structured meshes

  * Triangular meshes (impossible with the deal.II library, but adding a 
    different FEM library could make it possible; MEPLS would still depend on 
    deal.II tensor classes to work)

  * Internal parallelization of individual simulation runs

### Implementation affecting only the FEM solver
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
