### Particle filter estimator for multibody dynamics ###
This directory contains a C++ library for multibody dynamics estimation with 
a particle filter. Only planar dynamics models have been tested. 
For the theory behinds this work, refer to : 
  * J.L. Blanco, J.L. Torres, A. Gimenez-Fernandez, "Multibody dynamic systems as Bayesian Networks: applications to robust state estimation of mechanisms", Multibody System Dynamics, vol. 34, no. 2, pp. 103-128, 2015.  ([Draft PDF](http://ingmec.ual.es/~jlblanco/papers/blanco2015mds_bayesian_networks_DRAFT.pdf), [PDF](http://www.sciencedirect.com/science/article/pii/S0094114X15000919),  [BibTeX](http://ingmec.ual.es/aigaion2/index.php/export/publication/312/bibtex))

[![ScreenShot](https://raw.githubusercontent.com/MBDS/mbde-particle-filter/master/mbde-pf-screenshot.jpg)](https://www.youtube.com/watch?v=7Zru0oiz36g)

Compiling 
------------
Prerequisites: 
  * A C++ compiler
  * CMake
  * SuiteSparse
  * [MRPT](http://www.mrpt.org).

In Ubuntu, install them with: 

        sudo apt-get install build-essential cmake libsuitesparse-dev libmrpt-dev

You may need to install a newer version of MRPT [from here](http://www.mrpt.org/MRPT_in_GNU/Linux_repositories)

To build: 

        mkdir build
        cd build 
        cmake ..  
        # If you see no errors, go on, otherwise fix them!
        make        # To compile the library and examples
        make test   # To run unit tests

You can also build under Windows and Visual Studio.

Executing
------------
These programs come ready to be launched, and can be found in the `bin` 
directory after compiling:

  * `pf_test1`: One of the particle filter estimation experiments showed in the paper.
  * `ex_four_bars`: An example of a dynamic simulation of a four bar linkage. 

