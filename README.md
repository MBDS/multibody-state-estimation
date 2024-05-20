[![CI Linux](https://github.com/MBDS/multibody-state-estimation/actions/workflows/linux-build.yml/badge.svg)](https://github.com/MBDS/multibody-state-estimation/actions/workflows/linux-build.yml)
[![CI Check clang-format](https://github.com/MBDS/multibody-state-estimation/actions/workflows/check-clang-format.yml/badge.svg)](https://github.com/MBDS/multibody-state-estimation/actions/workflows/check-clang-format.yml)
[![Docs](https://readthedocs.org/projects/libmbse/badge/)](https://libmbse.readthedocs.io/)
[![Gitpod ready-to-code](https://img.shields.io/badge/Gitpod-ready--to--code-blue?logo=gitpod)](https://gitpod.io/#https://github.com/MBDS/multibody-state-estimation)

# The MultiBody State Estimation (MBSE) C++ library
A library for multibody forward and inverse dynamics, state estimation, and
parameter identification.

  * [C++ API docs](https://libmbse.readthedocs.io/).
  * [GitHub repo](https://github.com/MBDS/multibody-state-estimation).
  * [Applications](https://libmbse.readthedocs.io/en/latest/pageApps.html).

Related papers and approaches implemented in this repository:

## Paradigm #1: Based on factor-graphs
Forward and inverse dynamics, closed kinematic chains constraints solved as
factors in a factor graphs, and more.

For the theory behinds this work, refer to :
  * J.L. Blanco, A. Leanza, G. Reina, "A general framework for modeling and dynamic simulation of multibody systems using factor graphs", 2021.  ([ArXiV](https://arxiv.org/abs/2101.02874))

![ScreenShot](https://raw.githubusercontent.com/MBDS/multibody-state-estimation/master/doc/source/_static/mbde-fg-screenshot.png)

## Paradigm #2: Particle filter estimator
A Bayesian filter has been implemented as a particle filter.
Only planar dynamics models have been tested.
For the theory behinds this work, refer to :
  * J.L. Blanco, J.L. Torres, A. Gimenez-Fernandez, "Multibody dynamic systems as Bayesian Networks: applications to robust state estimation of mechanisms", Multibody System Dynamics, vol. 34, no. 2, pp. 103-128, 2015.  ([Draft PDF](http://ingmec.ual.es/~jlblanco/papers/blanco2015mds_bayesian_networks_DRAFT.pdf), [PDF](http://dx.doi.org/10.1007/s11044-014-9440-9),  [BibTeX](http://ingmec.ual.es/aigaion2/index.php/export/publication/289/bibtex))

<center> [Open video](https://www.youtube.com/watch?v=7Zru0oiz36g) </center>

![ScreenShot](https://raw.githubusercontent.com/MBDS/multibody-state-estimation/master/doc/source/_static/mbde-pf-screenshot.jpg)

## My config for running this project
Dependencies:
  * A C++ compiler
  * CMake
  * SuiteSparse
  * [MRPT](https://github.com/mrpt/mrpt/) (>=2.0.0)
    * Ubuntu: Use [this PPA](https://launchpad.net/~joseluisblancoc/+archive/ubuntu/mrpt)
  * [GTSAM](https://github.com/borglab/gtsam) =4.2
    * Ubuntu: gtsam was built by using system's default Eigen v3.4 with command : `cmake -DGTSAM_USE_SYSTEM_EIGEN=ON ..`

In Ubuntu, install dependencies with:

        sudo apt-get install build-essential cmake libsuitesparse-dev 

To build:

        mkdir build
        cd build
        cmake ..  
        # If you see no errors, go on, otherwise fix them!
        make        # To compile the library and examples
        make test_legacy   # To run unit tests

You should also be able to compile this project under Windows and Visual Studio.

To Run the FourBar.YAML example do the following in terminal

        
        <path to build directory>/bin/mbse-dynamic-simulation --mechanism <path_to_yaml_model_definition> [other_optional_arguments]
        

## Using mbse as a library in a user program

In your CMake project, add:

        find_package(mbse REQUIRED)
        # then in your target:
        target_link_libraries(${PROJECT_NAME} mbse::mbse)

For CMake to find the library, in `cmake-gui` or `ccmake`, set the variable `mbse_DIR` to
`BUILD_DIR/cmake`, where `BUILD_DIR` is the compilation directory where you built MBSE.

## License

Distributed under the 3-clause BSD license.
