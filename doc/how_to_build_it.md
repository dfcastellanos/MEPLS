
@page HowToBuild
 
# How to build
 
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

##Linux  

To build deal.II, unpack its sources in some directory `/path/to/dealii/sources`. Then, to 
configure and build it in, e.g., `/path/to/dealii/build`, you can do:

```sh
mkdir /path/to/dealii/build
cd /path/to/dealii/build
cmake /path/to/dealii/sources -DDEAL_II_WITH_UMFPACK=on
make # use the argument -j<N> for a build using N parallel threads
```
@note Alternatively, you can do `make install` to install the library in you system, but is not 
mandatory.

If you have any doubts about the process of building deal.II, you can find a more extensive 
explanation on its documentation (see [deal.II installation](https://www.dealii.org/current/readme
.html#installation))

Now, we can compile some of the MEPLS models. First, unpack MEPLS in some 
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

@note Instead of adding the flags `-DDEAL_II_DIR` and `-DMEPLS_DIR`, you can also set the 
corresponding environment variables as `export DEAL_II_DIR=/path/to/dealii/build` and `export 
MEPLS_DIR=/path/to/MEPLS_dir`. Also, if you chose to install deal.II in your system, CMake should
find it automatically.

After calling make, an executable named `run_sim` should appear in `/path/to/model/build`.

Each model is built using its own `CMakeLists.txt`, so different models might accept different 
configuration flags. However, the flag `DEFINE_DEBUG` turns on MEPLS debug mode,

  * `-DDEFINE_DEBUG=[on/off]` defines whether to build the model using the debug
    mode. In debug mode, MEPLS performs intensive checks, which can reduce the 
    performance significantly. If the flag is not set, the default value is `off`.
    

##Windows
Building on Windows has not been tested yet. You can follow [deal.II]'s building
instructions for Windows, after which building a MEPLS model should be straightforward.

##macOS
Building on macOS has not been tested yet. You can follow [deal.II]'s building
instructions for macOS, after which building a MEPLS model should be straightforward.

[deal.II]: https://www.dealii.org/
