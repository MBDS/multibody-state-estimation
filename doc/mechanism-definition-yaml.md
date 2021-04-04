\page pageMechDefYaml Mechanism definition YAML files

\tableofcontents

MBSE allows users to define mechanisms by means of a human-friendly YAML file.

Example: [mechanisms/fourbars1-with-rel-angle.yaml](https://github.com/MBDS/multibody-state-estimation/blob/master/config/mechanisms/fourbars1-with-rel-angle.yaml)
\include mechanisms/fourbars1-with-rel-angle.yaml

\section secStructure 1. YAML file structure

Prototype file structure:

	%YAML 1.2
	---
	# Optional top comment section.
	#
	# Define points (mandatory)
	points:
	  # 0 (="A")
	  - {x: 0, y: 0, fixed: true}
	  # 1
	  # ...
	# Define bodies (mandatory)
	planar_bodies:
	  # 0
	  - points: [0, 1]
	    # ...
	# Relative coordinates (Optional section)
	relative_coordinates:
	  # 0
	  - type: RelativeAngleAbsoluteDOF


\note YAML requires correct indentation to be observed

Required sections are:
  * `points`: one entry per defined "point", both moving and fixed ones.
  * `planar_bodies`: each planar body is defined by means of exactly *two* points.
    Length, mass, inertia with respect to the first point (`I0`),
    and center of gravity (`cog`) must be also specified, possibly using equations
    that get evaluated (in that order) substituting these variables:
    * `index`: 0-based index of the body,
    * `length`: the length.
    * `mass`: The mass.
    * `I0`: The inertia.


\section secAPI 2. C++ API

Usage:

\code
// Load mechanism model:
const auto yamlData = mrpt::containers::yaml::FromFile(filename);
const mbse::CModelDefinition model = mbse::CModelDefinition::FromYAML(yamlData);
\endcode
