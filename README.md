[![MBSE](https://circleci.com/gh/MBDS/multibody-state-estimation.svg?style=svg)](https://circleci.com/gh/MBDS/multibody-state-estimation)

# The MultiBody State Estimation (MBSE) C++ library
A library for multibody forward and inverse dynamics, state estimation, and
parameter identification.

Related papers and approaches implemented in this repository:

## Paradigm #1: Based on factor-graphs
Forward and inverse dynamics, closed kinematic chains constraints solved as
factors in a factor graphs, and more.

Refer to papers:

For the theory behinds this work, refer to :
  * J.L. Blanco, A. Leanza, G. Reina, "A general framework for modeling and dynamic simulation of multibody systems using factor graphs", 2021.  ([ArXiV](http://xxx)

![ScreenShot](docs/mbde-fg-screenshot.png)

## Paradigm #2: Particle filter estimator
A Bayesian filter has been implemented as a particle filter.
Only planar dynamics models have been tested.
For the theory behinds this work, refer to :
  * J.L. Blanco, J.L. Torres, A. Gimenez-Fernandez, "Multibody dynamic systems as Bayesian Networks: applications to robust state estimation of mechanisms", Multibody System Dynamics, vol. 34, no. 2, pp. 103-128, 2015.  ([Draft PDF](http://ingmec.ual.es/~jlblanco/papers/blanco2015mds_bayesian_networks_DRAFT.pdf), [PDF](http://dx.doi.org/10.1007/s11044-014-9440-9),  [BibTeX](http://ingmec.ual.es/aigaion2/index.php/export/publication/289/bibtex))

[![ScreenShot](docs/mbde-pf-screenshot.jpg)](https://www.youtube.com/watch?v=7Zru0oiz36g)


## Compiling
Dependencies:
  * A C++ compiler
  * CMake
  * SuiteSparse
  * [MRPT](http://www.mrpt.org) (>=2.0.0)
    * Ubuntu: Use [this PPA](https://launchpad.net/~joseluisblancoc/+archive/ubuntu/mrpt)
  * [GTSAM]()
    * Ubuntu: Use [this PPA](https://launchpad.net/~joseluisblancoc/+archive/ubuntu/gtsam-develop)

In Ubuntu, install dependencies with:

        sudo apt-get install build-essential cmake libsuitesparse-dev libmrpt-dev libgtsam-dev

To build:

        mkdir build
        cd build
        cmake ..  
        # If you see no errors, go on, otherwise fix them!
        make        # To compile the library and examples
        make test_legacy   # To run unit tests

You should also be able to compile this project under Windows and Visual Studio.

## Executing
These programs come ready to be launched, and can be found in the `bin`
directory after compiling:

  * `pf_test1`: One of the particle filter estimation experiments showed in the paper.
  * `ex_four_bars`: An example of a dynamic simulation of a four bar linkage.

## Using mbse as a library in a user program

In your CMake project, add:

        find_package(mbse REQUIRED)
        # then in your target:
        target_link_libraries(${PROJECT_NAME} mbse::mbse)

For CMake to find the library, in `cmake-gui` or `ccmake`, set the variable `mbse_DIR` to
`BUILD_DIR/cmake`, where `BUILD_DIR` is the compilation directory where you built MBSE.
