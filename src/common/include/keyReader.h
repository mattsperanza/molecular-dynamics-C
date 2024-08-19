// Author(s): Matthew Speranza
#ifndef KEYREADER_H
#define KEYREADER_H

#include "../system/system.h"

/**
 * Molecular Dynamics C Keywords:
 * - [] means optional (must complete all optional arguments if started)
 * - last definition of the same variable overwrites what is currently there
 * - all keywords are case-insensitive
 * - if only one axis is given all axis takes that length (in that case it must only contain one float)
 * - if only two axis are given the program will fail
 * - if only one boxDim is given, all dimensions will take on that length
 * - if only one grid count is given, all dimensions get that value
 * - grid counts must be powers of 2
 * - if grid counts are not given, they are automatically set
 * - all keywords have restrictions due to datatypes listed and are computer arch. dependent
 *
 *
 * verbose (bool) - logging level on or off
 * dt (float) - sets timestep (nanoseconds) (default 1) - overwrite potential
 * dtNano (float) - same as above - overwrite potential
 * dtAtto (long) - sets dynamics timestep (attoseconds) (default 1e3) - overwrite potential
 * printThermoEvery (long) - prints thermodynamic information every ? steps (default 1e4)
 * printRestartEvery (long) - prints restart files ever ? steps (default 1e4)
 * steps (long) - number of dynamics steps to take (default 1e9)
 * temp (float) - temperature of the system (kelvin) (default 298K)- overwrite potential
 * temperature (float) - same as above - overwrite potential
 * cutoff (float) - neighborlist cutoff (angstrom)
 * buffer (float) - neighborlist buffer (angstrom)
 * A-axis (float,[float,float]) - A-axis (Ax,[Ay,Az]) (angstrom) - overwrite potential
 * B-axis (float,[float,float]) - B-axis (Bx,[By,Bz]) (angstrom) - overwrite potential
 * C-axis (float,[float,float]) - C-axis (Cx,[Cy,Cz]) (angstrom) - overwrite potential
 * boxDim (float,[flaot,float]) - Cubic dimensions (X,Y,Z) (angstrom) - overwrite potential
 * pmeAlpha (float) - ewald paremeter (default ___)
 * pmeBeta (float) - ewald paremeter (default ___)
 * pmeOrder (int) - b-spline order (default 5)
 * pmeGridCount (int,[int,int]) - PME grid nodes (nX,nY,nZ)
 * polarization (char*) - Polarization type (mutual,direct,none)
 * forcefield (filepath) - path to force field file - overwrite potential
 * parameters (filepath) - same as above - overwrite potential
 * params (filepath) - same as above - overwrite potential
 * patch (filepath,[filepath,...]) - read in molecule parameters - overwrite potential if atomType number overlaps forcefield
 *
 */

static char* MD_C_Keywords[25] =
 {"verbose",
 "dt", "dtNano", "dtAtto",
 "steps",
 "printThermoEvery",
 "printRestartEvery",
 "temp", "temperature",
 "cutoff",
 "buffer",
 "a-axis","b-axis","c-axis","boxDim",
 "pmeAlpha",
 "pmeBeta",
 "pmeOrder",
 "pmeGridCount",
 "polerization",
 "forcefield", "parameters", "params",
 "patch",
 "printArchiveEvery"
};

void readKeyFile(System* system, char* keyFile);
void handleArgs(Vector* args, System* system);

#endif //KEYREADER_H
