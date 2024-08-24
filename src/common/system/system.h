// Author(s): Matthew Speranza
#ifndef SYSTEM_H
#define SYSTEM_H

#include <forceFieldReader.h>
#include <stdbool.h>
#include <vector.h>

#include "defines.h"

/**
 * Defines all the general system variables, computational details, and other information.
 * Pointer to the system struct is passed essentially everywhere.
 */
enum Polarization {NONE, DIRECT, MUTUAL};
typedef struct System {
 // Molecular System
 int nAtoms;
 REAL volume; // Angstroms^3
 REAL pressure; // atmospheres
 REAL density; // amu / ANG^3
 REAL particleDensity; // nAtoms / ANG^3
 REAL temperature; // Kelvin
 REAL** multipoles; // Force field definitions of multipolar charge distribution [nAtoms][cartesian multipole d.o.f. - 10 for now]
 REAL** frameDef; // Frame definition atoms for multipoles [nAtoms][5] - up to 4 atoms then an enum at the end for frame type
 REAL** rotatedMpoles; // Rotated multipoles
 int* atomTypes; // Atom forcefield type
 int* atomClasses; // Atom forcefield class
 int** polExclusions; // Exclusions for polarization [nAtoms][5] - up to 5 atoms
 Vector* verletList; // Indices in X of atoms within cutoff+buffer distance
 Vector* list12; // Indices in X of atoms every atom is bonded to vector of ints --> 1-2 lists
 Vector* list13; // Indices in X of atoms every atom is 1-3 bonded to vector of ints
 Vector* list14; // Indices in X of atoms every atom is 1-4 bonded to vector of ints
 Vector* list15; // Indices in X of atoms every atom is 1-5 bonded

 Vector bonds; // All 2-len paths in structure (bonds) (vector of int[2])
 Vector piTorsions; // Indices in bonds array for pi-torsions
 Vector angles; // All 3-len paths in bonds (angles) (vector of int[3])
 Vector strBendAngleIndex; // Index into angles array for stretch-bend angles
 Vector urayBradIndex; // Index into angles array for Urey-Bradley angles
 Vector oopBendIndex; // Index into angles array for out-of-plane bending angles
 Vector torsions; // All 4-len paths in structure in bond order (torsions) (vector of int[4])
 Vector torsionTorsion; // All 5-len paths in structure in bond order (tortors) (vector of int[5])

 REAL boxDim[3][3]; // Box axis definitions (ATM) [A,B,C][x,y,z]
 REAL minDim[3]; // Minimum box dimensions (ANG) [x,y,z]
 char** atomNames; // Atom periodic table name [nAtoms][name]
 REAL* protons; // Number of protons [nAtoms]
 REAL* valence; // Number of valence electrons [nAtoms]

 // Simulation method parameters
 int dtAtto;
 REAL dtFemto;
 long steps; // Number of dynamics steps to run
 long currentStep; // Current simulation time in attoseconds
 long printThermoEvery; // Print energy information every n femto-seconds
 long printRestartEvery; // Print restart *.dyn file every n femto-seconds
 long printArchiveEvery; // Print snap into *.arc file every n femto-seconds
 REAL ewaldAlpha; // Gaussian parameter
 REAL ewaldBeta; // Gaussian parameter
 REAL ewaldOrder; // Order of b-splines
 int* pmeGridspace; // number of gridpoints in each dimension [nX, nY, nZ]
 REAL*** pmeGrid; // Grid of splined multipoles [nX][nY][nZ]
 REAL* pmeGridFlat; // Grid of splined multipoles [nX*nY*nZ]
 REAL vdwCutoff; // Neighborlist cutoff in angstroms
 REAL vdwTaper; // Tapering percentage of vdw cutoff
 REAL ewaldCutoff; // Ewald cutoff in angstroms
 REAL realspaceBuffer; // Addtion to cutoff to buffer neighborlist builds
 ForceField* forceField; // Force field definitions
 enum Polarization polarization; // Polarization for amoeba

 // Integeration Variables
 int nDOF; // nAtoms + nActiveLambdas;
 REAL* DOF; // Degrees of freedom for integration [X, activeLambdas]
 REAL* DOFFrc; // Force on degrees of freedom
 REAL* X; // Interleaved atomic positions (ANG) (x,y,z) [nAtoms*3]
 REAL* M; // Atomic masses [nAtoms]
 REAL* V; // Interleaved atomic velocity (ANG/ns) (Vx, Vy, Vz) [nAtoms*3]
 REAL* A; // Interleaved atomic accelerations (ANG/ns^2) (Ax, Ay, Az) [nAtoms*3]
 REAL* F; // Interleaved atomic forces (Fx,Fy,Fz) (kcal/mol/ANG) [nThread][nAtoms*3]
 REAL* lambdas; // Atom lambdas between 0-1 [nAtoms]
 REAL* thetas; // Atom thetas - converted into lambda 0-1 [nLambdaVariables]
 REAL* thetaM; // Theta masses
 REAL* thetaV; // Atom theta velocities
 REAL* thetaA; // Atom theta acceleration
 REAL* thetaF; // Force on theta
 int* activeLambdas; // which atom index have important lambda values (-1 if all)
 int nActiveLambdas; // Length of previous array

 // Computer definitions
 bool verbose;
 char* structureFileName; // Can't deallocate since user input
 bool isPDB;
 bool isXYZ;
 char* remark; // First line of xyz file that contains the atomnumber
 char* structureFilePath; // where all output file writing is directed and restart files should be located
 char* forceFieldFile; // Path to force field
 Vector patchFiles; // Vector of char* indicating patch files
 char* keyFileName; // can be located anywhere -> useful to set up script one time and execute many // Cant deallocate
 int nThreads; // = 1; // number of threads assigned to this system
 int* threadIDs[1]; // the new id's assigned to threads of this system
} System;

/**
 * Applies PBC.
 *
 */
REAL pbc(REAL dx, REAL axisLen);


#endif //SYSTEM_H
