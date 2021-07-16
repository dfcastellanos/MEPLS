

@page mainpage


<div id="TextBox">

## About 

MEPLS (Mesoscale Elasto-Plasticity Simulator) is an open-source C++ framework 
for simulating the plastic deformation of materials based on the combination of 
stochastic processes and solid mechanics. Its goal is to provide modular, 
efficient, and easy-to-use tools to build simulations for a broad range of
physical scenarios. MEPLS abstracts away unnecessary complexity associated with the model 
implementation, allowing the user to focus on the physics of the problem. To achieve these goals,
MEPLS is developed around three pillars:

  * **Efficient** — it is fully implemented in C++, with a back-end based on 
    well-established efficient libraries such as deal.II and UMFPACK.

  * **User-friendly** — it is easy to build and has minimal dependencies. It 
    is well-documented and has a gallery of examples.

  * **Extensible** — it is object-oriented and template-based. Its modularity 
    and simplicity facilitate the task of extending it according to the user's 
    needs.

MEPLS enables the user to study the impact that microstructural properties have on the macroscale
mechanical behavior of a material. This question is especially challenging to answer when 
those properties are spatially heterogeneous and statistically distributed, which is, in turn, 
the most common situation (see [Background] (@ref Background)). MEPLS allows to study 
the influence that many different model features or ingredients have on a material's mechanical 
properties. Example of such features are:

* the statistics of local elastic properties, local yield criteria, plastic deformation, structural 
evolution, etc.
* the driving protocol: strain-controlled tests, creep, cyclic loading, etc.
* the loading mode: simple shear, compression, bending, etc.

You can find an extensive list in the [features section](@ref Features).


## Tutorial 

Each [tutorial](@ref Tutorial) step introduces different tools and illustrates how to use them. 
After going through the tutorial steps, you will have a good idea of how to use MEPLS for implementing you 
ouwn model.

## Gallery

New models can be implemented following the tutorial and the MEPLS documentation. However, 
sometimes a model might require a use of MEPLS more advanced than illustrated in the tutorial. In
this case, the [gallery](@ref Gallery) of models shows advanced examples built for specific 
purposes. 


## Documentation
This website provides you with an interactive MEPLS [documentation](@ref Documentation). The 
**search box** on the top right corner allows you to search through the documentation for 
functions, classes, etc. The source code in this website is also interactive, allowing you to
click on it and jump to its documentation.


## Download

You can download the source files from the [project's repository] on GitHub.

## Contact
MEPLS is developed and maintained by David Fernández Castellanos. You can report issues and bugs 
in the [project's repository](https://github.com/kastellane/MEPLS). You can contact the author 
through the methods provided on the [author's website] for longer discussions regarding, e.g., 
requests and ideas for the project or if you need some help to use MEPLS.


## License
MEPLS is open source. You can freely use it, redistribute it, and/or modify it
under the terms of the [Creative Commons Attribution 4.0 International Public 
License](https://creativecommons.org/licenses/by/4.0/).


##Last update
Jul 10 2021


<br></div> 


[deal.II]: https://www.dealii.org/
[project's website]: https://mepls.davidfcastellanos.com/
[author's website]: https://www.davidfcastellanos.com/contact
[project's repository]: https://github.com/kastellane/MEPLS

