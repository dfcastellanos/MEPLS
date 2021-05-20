

@page mainpage

### About

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


MEPLS provides easy-to-use tools to build mesoscale models of the plastic 
deformation of materials and to statistically sample them with Markov chain 
Monte Carlo methods. MEPLS enables the user to study the impact that 
microstructural properties have on the macroscale elasto-plastic behavior of a 
material. This question is especially challenging to answer when those properties
are spatially heterogeneous and statistically distributed, which is, in turn, the
most common situation. The built-in MEPLS tools abstract away unnecessary complexity 
associated with model specification, allowing the user to focus on the physics 
of the problem. Thus, the user can easily study the influence of: the local 
statistics of elastic constants, yield criteria, plastic deformation, structural
evolution, etc.; the driving protocol: strain-controlled tests, creep, cyclic 
loading, etc.; the loading mode: simple shear, compression, bending, etc. 
Moreover, suppose a feature is not readily available. In that case, chances are
that the user can seamlessly implement it thanks to MEPLS extensibility, 
availability of examples, and thorough documentation.

MEPLS allow to study the influcence of many different ingredients such as:
* local statistics of elastic constants, yield criteria, plastic deformation, structural evolution, etc.
* driving protocol: strain-controlled tests, creep, cyclic loading, etc.
* loading mode: simple shear, compression, bending, etc.

You can find an extensive list in the features section.


### Download

You can download the source files from the [project's repository] on GitHub.


### Contact

You can report issues and bugs in the [project's repository]. You can contact the author through 
the methods provided on the [author's website] for longer discussions regarding, e.g., requests 
and ideas for the project or if you need some help to use MEPLS.




[deal.II]: https://www.dealii.org/
[project's website]: https://mepls.davidfcastellanos.com/
[author's website]: https://www.davidfcastellanos.com/
[project's repository]: https://github.com/kastellane/MEPLS

