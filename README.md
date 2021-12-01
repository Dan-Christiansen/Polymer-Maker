# Polymer-Maker
Package to create linear zwitterionic homopolymer structure and topology files for use with GROMACS OPLS-AA molecular dynamics simulations.
### Usage
---
Import _polymer_maker_ with
> from polymer_maker import polymer_maker

Create a polymer structure and system topology with
> polymer_maker(n_polymers, angle)

_Inputs_

* polymer_formula - Text file describing the polymer to be made. See sample file for details
* n_polymers - Number of polymers to be in the final system
* angle - Degrees of rotation between neighboring monomers

_Outputs_

* newmol.gro - Structure file of polymer
* topol.top - Topology file for system
* newmol.itp - Topology file of polymer
* posre.itp - Atom position restraint file. Unimportant for most applications

### Known issues:
---
- Currently dependent on a fixed path to _gmx_ (search for /usr/local/gromacs_gmx/bin/gmx in code)

Made by Daniel Christiansen
For Molecular Simulations Laboratory, University of Illinois at Chicago, Chicago, IL 60608
