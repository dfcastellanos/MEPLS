
@page HowToBuild
 
####Linux
MEPLS uses CMake for building and depends on the [deal.II] library for using the
Finite Element Method. Compatible deal.II versions are between the 8.5 
and the 9.1. The recommended is the 9.0, which can be found in 
[this repository](https://github.com/dealii/dealii/tree/dealii-9.0). We will need
the LAPACK library, which is very standard, and you might already have it 
installed in your system. If not, you can install it from your OS distribution 
repositories. (For example, for Ubuntu-based distributions, simply do
`sudo apt-get install liblapack-dev`).

To build deal.II, unpack its sources in some directory `/path/to/dealii/sources`. 
Then, to configure and build it in, e.g., `/path/to/dealii/build`, do:

```sh
mkdir /path/to/dealii/build
cd /path/to/dealii/build
cmake /path/to/dealii/sources -DDEAL_II_WITH_UMFPACK=on
make # use the argument -j<N> for a build using N parallel threads
make test
```

Now, we can compile some of the MEPLS models. First, unpack MEPLS in some 
directory `/path/to/MEPLS_dir`. Let’s build the model located in 
`/path/to/MEPLS_dir/models/example`. To build it in `/path/to/model/build`, simply
do:

```sh
mkdir /path/to/model/build
cd /path/to/model/build
cmake /path/to/MEPLS_dir/models/example -DDEAL_II_DIR=/path/to/dealii/build
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
instructions for Windows, after which building a MEPLS model should be straightforward.

####macOS
Building on macOS has not been tested yet. You can follow [deal.II]'s own building
instructions for macOS, after which building a MEPLS model should be straightforward.

[deal.II]: https://www.dealii.org/
