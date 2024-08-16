// Author(s): Matthew Speranza
#ifndef SYSTEM_H
#define SYSTEM_H

#include <stdbool.h>
#include <stdlib.h>
#include <vector.h>

#include "defines.h"

/**
 * Defines all of the general system variables, computational details, and other information.
 * Pointer to the system struct is passed essentially everywhere.
 */

enum ForceField {AMOEBA, CHARMM};
typedef struct System {
 // Molecular System
 int nAtoms;
 REAL volume; // Angstroms^3
 REAL pressure; // atmospheres
 REAL density; // amu / ANG^3
 REAL temperature; // Kelvin
 REAL** multipoles; // Force field/Quantum definitions of multipolar charge distribution [nAtoms][cartesian multipole d.o.f.]
 int* atomTypes; // Atom forcefield type
 Vector** bondList; // Indices in X of atoms every atom is bonded to vector of ints
 REAL** boxDim; // Box axis definitions (ATM) [A,B,C][x,y,z]
 char** atomNames; // Atom periodic table name [nAtoms][name]
 REAL* protons; // Number of protons [nAtoms]
 REAL* valence; // Number of valence electrons [nAtoms]

 // Simulation method parameters
 int dtAtto;
 REAL dtFemto;
 long currentStep; // Current simulation time in attoseconds
 REAL ewaldAlpha; // Gaussian parameter
 REAL ewaldBeta; // Gaussian parameter
 REAL ewaldOrder; // Order of b-splines
 REAL* pmeGridspace; // number of gridpoints in each dimension [nX, nY, nZ]
 REAL*** pmeGrid; // Grid of splined multipoles [nX][nY][nZ]
 REAL* pmeGridFlat; // Grid of splined multipoles [nX*nY*nZ]
 REAL realspaceCutoff; // Neighborlist cutoff in angstroms
 REAL realspaceBuffer; // Addtion to cutoff to buffer neighborlist builds
 enum ForceField ff; // Which force field is being used

 // Integeration Variables
 int nDOF; // nAtoms + nActiveLambdas;
 REAL* DOF; // Degrees of freedom for integration [X, activeLambdas]
 REAL* DOFFrc; // Force on degrees of freedom
 REAL* X; // Interleaved atomic positions (ANG) (x,y,z) [nAtoms*3]
 REAL* F; // Interleaved atomic forces (Fx,Fy,Fz) (kcal/mol/ANG) [nThread][nAtoms*3]
 REAL* lambdas; // Atom lambdas between 0-1 [nAtoms]
 REAL* thetas; // Atom thetas - converted to lambda 0-1 [nLambdaVariables]
 REAL* thetaFrcs; // Force on theta
 int* activeLambdas; // which atom index have important lambda values (-1 if all)
 int nActiveLambdas; // Length of previous array

 // Computer definitions
 char* structureFileName;
 bool isPDB;
 bool isXYZ;
 char* remark; // First line of xyz file that contains the atomnumber
 char* structureFilePath; // where all output file writing is directed and restart files should be located
 char* keyFileName; // can be located anywhere -> useful to set up script one time and execute many
 int nThreads; // = 1; // number of threads assigned to this system
 int* threadIDs[1]; // the new id's assigned to threads of this system
} System;

#endif //SYSTEM_H
