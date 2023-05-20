
@page Background

<div id="TextBox">

# Background


## Motivation
The history of innovation is linked to the discovery of novel materials with superior properties. From wood to rock, and from rock, to concrete and steel. From glass and metals to plastics. Nowadays, the production of a superior material usually results from extended and combined research endeavors. Regarding mechanical properties --that is, a material’s response to applied stresses-- materials scientists and physicists have looked into ever lower length scales to understand their microscopic origins. Mechanical properties are the consequence of the so-called microstructure: the material’s structure at a scale small enough such that its relevant constituents (be it atoms, molecules, grains, fibers, etc.) can be perceived. At the other end of the spectrum, at the macroscale, we perceive a material's response to stress as elastic, plastic, viscoelastic etc. These macroscale concepts have been described, combined, and exploited since the late 19th century. Nowadays, a primary goal is to establish links between microstructural properties and macroscale behavior. Such links can guide the design of new materials that meet the demands of specific applications and industries. On the other hand, they can enhance our understanding of the small-scale physical mechanisms ultimately responsible for the macroscale, observable behavior.


## Why mesoscale?
Historically, plastic deformation has been described as a smooth and deterministic process akin to the laminar flow of a fluid. However, modern technology has enabled researchers to systematically investigate materials’ behavior below the micron scale. A widely different picture becomes apparent at such small scales: plastic deformation is highly non-linear, intermittent, exhibits stochastic features, and gives rise to spatial patterns. Qualitatively similar behaviors have been found in a wide diversity of materials such as, e.g., crystals, glasses, rocks, wood, paper, solid foams, or some complex fluids, among others. The fact that such a diversity of materials share similarities suggests that their plastic deformation results from laws operating at a scale at which microscopic differences between them are not essential. In other words, at a range of scales between the micro and the macroscale, a material’s plastic deformation can be described with great independence of microscopic details. If the scale is not too large, such that we can still perceive spatial variability, that scale is known as the mesocale.

Models operating at the mesoscale are an excellent tool for establishing the mentioned links between micro and macroscale descriptions of plastic deformation. Moreover, from a computational performance point of view, they can reach larger spatial and temporal scales than microscale methods. At the same time, they can reproduce fluctuating phenomena missing in macro-scale approaches. Such fluctuations, which can be misinterpreted as noise, are of the utmost importance since their statistics encode most of the information regarding the nature of the plastic deformation process.

<br>
<br>
<center>
<img alt="Logo" src="scales.png" width="100%"/>
</center>
<br>


## Designing a model
The guiding principle of mesoscale modeling is to replace the quasi-infinite details present at the microscale for a coarser, more manageable description. The coarser description must retain only the fundamental features necessary for describing the process of interest, in this case, plastic deformation. To reduce the complexity of the description, mesoscale models divide a material sample into mesoscale subdomains, known as elements. The elements have only a few internal variables, defined in terms of continuum mechanics. This is to be compared with, e.g., models based on the position and velocities of the tens or hundreds of atoms that would inhabit a single mesoscale element. The behavior of the elements is defined according to a set of simple local rules, representing the most relevant characteristics of the material's microstructure and plastic deformation mechanisms. By doing so, the mesoscale description suffers a loss of information with respect to the microscale one. Such loss means that we have imperfect knowledge about the elements' state. To capture such statistical uncertainty, the rules governing the behavior of the elements correspond to stochastic processes.

Lastly, plastic deformation is represented as a series of individual unit deformation events, or quanta, localized both in space and time. The effects of an event are modeled by the plastic deformation of a single element. Since the elements compose a solid material, they must obey elementary rules of solid mechanics, such as stress equilibrium (i.e., the continuum version of Newton's 3rd law). Stress equilibrium implies that when an element deforms, it influences all other elements due to the material's elastic behavior. This elastic influence leads to a highly correlated system that exhibits complex behavior. Eventually, the self-organization of the elements gives rise to the emergence of, e.g., spatial patterns and well-known macroscale laws of plastic deformation.


## Interpretating mesoscale models
Mesomodels divide a material sample into discrete mesoscale subdomains or elements and implement the physics of the model as a set of simple local rules. There are two deeply connected interpretations of such models:

* The first interpretation is that of a cellular automaton. Or, more accurately, a stochastic one. Cellular automata are well-known sandbox models. With them, the user can easily modify state transition rules to see the impact on the self-organization of the system and the emergent global behavior.
	
* On the other hand, from a purely statistical perspective, these models represent a Markov chain
 Monte Carlo (MCMC) method applied to the plastic deformation of materials. The goal of the model
  is to help find the probability distributions that describe a particular model’s outcome. 
  However, finding such distributions is a non-trivial task since they depend on the interplay of
   many highly correlated variables. A way to estimate them is using MCMC methods. MCMC depends 
   on rules which define the transition from one state to the next. By starting with an initial 
   configuration and iterating the transition rules, we obtain a sequence of random samples that 
   correspond to the system’s state over time. These samples allow us to estimate the desired 
   distributions. In this case, the MCMC transition rules correspond to the local mesoscale rules
    mentioned above, which model physical mechanisms and govern the behavior of the elements.


## Why MEPLS?
MEPLS is a framework that provides easy-to-use tools to build mesoscale models of the plastic 
deformation of materials under a broad range of physical scenarios. Moreover, MEPLS allows 
sampling a mesoscale model’s configuration space with Markov chain Monte Carlo methods. In this 
way, the user can study the impact that microstructural properties have on the macroscale 
mechanical behavior of a material. This question is especially challenging to answer when 
those properties are spatially heterogeneous and statistically distributed, which is, in turn, 
the most common situation. MEPLS provides a set of tools that abstract away unnecessary complexity 
associated with the model implementation, allowing the user to focus on the physics of the problem.

## References
@cite DFCastellanos_CRP @cite BudrikisNatCom @cite Talamali2012 @cite nicolas_deformation_2018 
@cite karimi_role_2016 @cite bulatov_stochastic_1994-1

<br></div>
