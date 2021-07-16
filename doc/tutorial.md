

@page Tutorial

<div id="TextBox">

# Tutorial

The following tutorial steps introduce MEPLS, after which you will have a good idea
of where to start to implement your model. Each step builds on top of the previous one, and the 
complexity gradually increases by adding new features. Sometimes, however, there is overlap 
between steps, since instead of adding a new feature, a previous one is re-implemented in a 
better way. Each tutorial step  has the following structure:

* **Introduction**: a brief description of what the tutorial step is about
* **The commented program**: the program, split into short pieces of code with detailed descriptions
* **Results**: the results of the model, how to visualized them, and how to interpret them
* **The complete program**: the step's source code
    
For the sake of clarity, the tutorial guides you through a specific use case. Specifically, the 
process of implementing an elasto-plastic mesoscale model of a material undergoing plastic 
deformation under a slow applied strain rate. This model aims to allow the user to tune 
microstructural properties (at the mesoscale) and study the impact on the stress-strain response 
and strain localization, such as brittle-to-ductile transition or shear-banding. The steps cover
 the following topics: 

* @ref Step1 
    * Introduction to the element and slip classes
    
* @ref Step2 
    * The elasticity solver
    * The system class
    * Plastic and driving events
    
* @ref Step3 
    * Modifying the external load
    * Triggering slip events based on rules
    * Stopping criterion
    * Evolution history
    
* @ref Step4 
    * The same dynamics as in Step 3, but using MEPLS built-in tools: 
        * Extremal dynamics 
        * Mechanical relaxation
    * Writing output data    
    * Animating the results
    
* @ref Step5 
    * Thermal slip events and the Kinetic Monte Carlo method
    * Combining different driving protocols
    * Running several repetitions in parallel
    
* Step 6 (to be added):
    * How to implement new elements and slips    

The tutorial source code is located in `MEPLS/tutorial`, and each step can be built like other 
MEPLS models, as explained in [How to build](@ref HowToBuild).            
     
@note After the tutorial, you can take a look at the @ref Gallery of models. Such models are
less documented, but they are natural extensions of the tutorial steps, illustrating how to 
implement particular model features. 

<br></div>




