// Author(s): Matthew Speranza
#ifndef DEFINES_H
#define DEFINES_H

// Try not to include anything in this file
typedef double REAL;
#define BOHR_TO_ANGSTROM 0.52917721090 // From FFX
#define ANGSTROM_TO_BOHR (1.0 / BOHR_TO_ANGSTROM)
#define BOHR_TO_ANGSTROM2 (BOHR_TO_ANGSTROM * BOHR_TO_ANGSTROM)
#define ANGSTROM_TO_BOHR2  (1.0 / BOHR_TO_ANGSTROM2)
// Coulomb Constant in units: kcal*Ang/(mol*electron^2)
#define ELECTRIC 332.0637133
#define HARTREE_TO_KCAL_PER_MOL = 627.5094738898777
#endif //DEFINES_H