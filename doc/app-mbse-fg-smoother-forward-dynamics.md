\page pageApp-mbse-fg-smoother-forward-dynamics mbse-fg-smoother-forward-dynamics

The program `mbse-fg-smoother-forward-dynamics` is a demo for forward
dynamic simulation using factor graphs (FG).

\section sec1 Examples of use

	mbse-fg-smoother-forward-dynamics \
	  --mechanism ../config/mechanisms/fourbars1.yaml \
	  --dt 5e-3 --end-time 2.0 --lag-time 15e-3
	
	mbse-viewer \
	  --mechanism ../config/mechanisms/fourbars1.yaml \
	  -q q.txt

See also: \ref pageMechDefYaml

\section sec2 CLI options

USAGE:
\verbatim

   mbse-fg-smoother-forward-dynamics  [-v] [--final-batch]
                                        [--show-factor-errors]
                                        [--dont-add-dq-constraints]
                                        [--dont-add-q-constraints]
                                        [--output-prefix <prefix>]
                                        [--smoother-iterations <xxx>]
                                        [--lag-time <xxx>] [--end-time
                                        <xxx>] [--dt <xxx>] --mechanism
                                        <YAML model definition> [--]
                                        [--version] [-h]


Where:

   -v,  --verbose
     Verbose console output

   --final-batch
     Run an additional final batch optimizer

   --show-factor-errors
     Show factor errors for the final state

   --dont-add-dq-constraints
     Do NOT add the dq manifold constraint factors

   --dont-add-q-constraints
     Do NOT add the q manifold constraint factors

   --output-prefix <prefix>
     Output files prefix

   --smoother-iterations <xxx>
     Smoother optimization maximum iterations

   --lag-time <xxx>
     Smoother lag time

   --end-time <xxx>
     Simulation end time

   --dt <xxx>
     Step time

   --mechanism <YAML model definition>
     (required)  Mechanism model YAML file

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.

\endverbatim
