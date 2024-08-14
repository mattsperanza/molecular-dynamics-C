# Script Directory
This directory contains mostly small methods for reading in user input.

## Sub-Directories
None

## Files
### md.c
Program entry point with main() function. Times execution time of the command.
### commandInterpreter.c
Parses user input and calls the appropriate function. Parses files and associates structure with force field. Prints logo.
### energy.c
Makes calls to create neighbor lists and calculate energies/forces.
### dynamics.c
Makes calls to create neighbor lists, calculate energies/forces, and integrate the equations of motion.

## Notes
- How should we differentiate between classical and quantum simulations?
    Through input alone. The energy.c should have an energyClassical method (also used by dynamics)
    and an energyQuantum method (called externally through a common energy method). The specific 
    method that is called can be determined once at the beginning.