# Parser Directory
This directory contains structure file parsing methods that fill out a system state.

## Sub-Directories
None

## Files
#### Each parser is required to fill out the system-state as fully as possible and return a pointer to the system.
#### Each parser is required to provide methods that deallocates all the memory that it allocates.
### key.c
Parses .key files that read "molecular-dynamics-c" keywords for simulation parameters.
### pdb.c
Parses .pdb files.
### xyz.c
Parses .xyz files.