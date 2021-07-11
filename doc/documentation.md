

@page Documentation

<div id="TextBox">


# Documentation

In the *Documentation > Namespace List* and *Documentation > Class List* tabs, you can access the
interactive MEPLS documentation. Also, the **search box** on the top right corner of this  
website allows you to search through the MEPLS documentation for namespaces, functions, classes 
or any other entity.

<br>

When navigating through the [tutorial](@ref Tutorial) steps of the [gallery](@ref Gallery) models, 
the programs' source code is also interactive, allowing you to click on different parts of the 
code and jump to the its documentation.


## Structure
The core of MEPLS is the header-only library located in `/path/to/MEPLS_dir/inc/MEPLS`. The MEPLS
tools are all defined in the namespace `mepls`. Within that namespace, there are several 
namespaces grouping classes and functions with similar responsibilities:


<br>
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
| `utils` | contains utilities that come in handy in different situations |

</center>

<br>
</div>