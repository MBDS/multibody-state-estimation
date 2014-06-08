### Example of particle filter for multibody dynamics ###

This directory contains a C++ library for multibody dynamics estimation with 
a particle filter. Only planar dynamics models have been tested.

[![ScreenShot](https://raw.github.com/jlblanco/mbde/master/mbde-pf-screenshot.jpg)](https://www.youtube.com/watch?v=7Zru0oiz36g)

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

