\page pageApp-mbse-dynamic-simulation mbse-dynamic-simulation

The program `mbse-dynamic-simulation` is a demo for forward dynamic simulation
solvers.

\section sec1 Examples of use

    bin/mbse-dynamic-simulation --mechanism ../config/mechanisms/fourbars1.yaml

On the GUI, press `+` or `-` to increase or reduce the simulation speed.
Close the GUI to end the program and see simulation statistics.

See also: \ref pageMechDefYaml

\section sec2 CLI options

\verbatim

USAGE:

   bin/mbse-dynamic-simulation  --mechanism <YAML model definition> [--]
                                [--version] [-h]


Where:

   --mechanism <YAML model definition>
     (required)  Mechanism model YAML file

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.


   mbse-dynamic-simulation

\endverbatim
