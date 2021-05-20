

@page Structure

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

<center>

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

</center>


#### Existing models & key features
The bundled models are located in `/path/to/MEPLS_dir/models/`. The available 
models are:

  * ESPCI:

    * Random local slip planes as glass microstructure
    
    * Monte Carlo simulation of a super-cooled liquid
    
    * Athermal quasistatic strain-controlled test
    
    * Unloading and reloading
    
    * Sampling local elasto-plastic response using the patch method
    



