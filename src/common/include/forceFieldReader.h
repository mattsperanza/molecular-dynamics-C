// Author(s): Matthew Speranza
#ifndef FORCEFIELDREADER_H
#define FORCEFIELDREADER_H
#include <vector.h>
#include "../system/defines.h"

/**
 * This was adapted from ForceFieldX's "ForceFieldFilter.java" class.
 */

enum ForceFieldName {AMOEBA, CHARMM};
enum ForceFieldParams {
  ATOM, ANGLE, ANGLEP, ANGTORS, BIOTYPE,
  BOND, CHARGE, MULTIPOLE, OPBEND, STRBND,
  PITORS, IMPTORS, STRTORS, TORSION, IMPROPER,
  TORTORS, UREYBRAD, VDW, VDW14, VDWPR, VDWPAIR,
  POLARIZE, RELATIVESOLV, SOLUTE
};
static char* ForceFieldParamsStr[24] = {
  "ATOM", "ANGLE", "ANGLEP", "ANGTORS", "BIOTYPE",
  "BOND", "CHARGE", "MULTIPOLE", "OPBEND", "STRBND",
  "PITORS", "IMPTORS", "STRTORS", "TORSION", "IMPROPER",
  "TORTORS", "UREYBRAD", "VDW", "VDW14", "VDWPR", "VDWPAIR",
  "POLARIZE", "RELATIVESOLV", "SOLUTE"
};
typedef struct Atom {
  int type, aClass, atomicNum, valence;
  REAL atomicMass;
  char name[4], environment[25]; // end with \0
} Atom;
enum AngleMode {NORMAL, IN_PLANE}; // Supported angle def.
enum AngleFunction {ANGLE_HARMONIC, SEXTIC}; // Supported angle functions
typedef struct Angle {
  int aClasses[3];
  REAL forceConstant;
  REAL angle[4];
  enum AngleMode angleMode;
  enum AngleFunction angleFunction;
} Angle;
typedef struct AngTors {
  int aClasses[4];
  REAL forceConstants[6];
} AngTors;
typedef struct BioType {
  int index, atomType;
  char atomName[5]; // end with \0
  char moleculeName[100]; // end with \0
} BioType;
enum BondFunction {BOND_HARMONIC, QUARTIC, FLAT_BOTTOM_HARMONIC, FLAT_BOTTOM_QUARTIC};
typedef struct Bond {
  int atomClasses[2];
  REAL forceConstant, distance, flatBottomRadius;
  enum BondFunction bondFunction;
} Bond;
// Includes charge as well
enum MultipoleFrameDef {MPOL_NONE, ZONLY, ZTHENX, BISECTOR, ZTHENBISECTOR, THREEFOLD};
typedef struct Multipole {
  REAL multipole[10]; // quadrupole * 1/3, off diag * 2/3
  int frameAtomTypes[4];
  enum MultipoleFrameDef frameDef;
} Multipole;
typedef struct OPBend {
  int atomClasses[4];
  REAL forceConstant;
} OPBend;
typedef struct StrBend {
  int atomClasses[3];
  REAL forceConstants[2];
} StrBend;
typedef struct PiTors {
  int atomClasses[2];
  REAL forceConstant;
} PiTors;
typedef struct ImpTors {
  int atomClasses[4];
  REAL forceConstant, phase, periodicity;
} ImpTors;
typedef struct StrTors {
  int atomClasses[4];
  REAL forceConstants[9];
} StrTors;
enum TorsionMode {TORS_NORMAL, TORS_IMPROPER};
typedef struct Torsion {
  int terms;
  int atomClasses[4], periodicity[3];
  REAL amplitude[3], phase[3];
  enum TorsionMode torsionMode;
} Torsion;
typedef struct TorTors {
  int atomClasses[5], gridPoints[2];
  REAL torsion1[625], torsion2[625], energy[625];
} TorTors;
typedef struct UReyBrad {
  int atomClasses[3];
  REAL forceConstant, distance;
} UReyBrad;
enum VdWType {VDW_NORMAL, VDW_14};
typedef struct VdW {
  int atomClass;
  REAL radius, wellDepth, reductionFactor;
  enum VdWType vdwType;
} VdW;
typedef struct VdWPair {
  int atomClasses[2];
  REAL radius, wellDepth;
} VdWPair;
// Currently isotropic, xx=yy=zz && mult with dirac delta
typedef struct Polarize {
  int atomType;
  int polarizationGroup[6];
  REAL polarizabilityTensor[3][3];
  REAL thole, ddp;
} Polarize;
typedef struct RelativeSolv {
  char residueName[10]; // end with \0
  REAL solvationEnergy;
} RelativeSolv;
typedef struct Solute {
  int atomType;
  REAL diameters[3]; // 0 = p-b, 1 = cos, 2 = gk
  REAL sneck;
} Solute;
// Defines all atom types and interactions between atom types
typedef struct ForceField {
  enum ForceFieldName name;
  Vector atom; // pointer to single vector with atom struct elements
  Vector angle;
  Vector angleParams; // contains angle params for every angle in System struct
  Vector angTors;
  Vector bioType;
  Vector bond;
  Vector bondParams; // Bond ff param definitions for every 12 path - vector of Bond*
  Vector multipole;
  Vector opBend;
  Vector opBendParams; // Out-of-plane bend ff param definitions for every opBend angle
  Vector strBend;
  Vector strBendParams; // Stretch-bend ff param definitions for every strBend angle
  Vector piTors;
  Vector piTorsParams; // Pi-torsion ff param definitions for every 12 path
  Vector impTors;
  Vector strTors;
  Vector torsion;
  Vector torsionParams; // Torsion ff param definitions for every 14 path
  Vector torTors;
  Vector tortorParams; // Tortor ff param definitions for every 15 path
  Vector uRayBrad;
  Vector urayBradParams;
  Vector vdw;
  Vector vdwPair;
  int nAtoms; // same as in system just for deleting
  REAL vdwN, vdwM, vdwDelta, vdwGamma;
  REAL* reductionFactors; // [nAtoms]
  REAL** wellDepths; // [nAtoms][neighborList.len]
  REAL** rMin; // [nAtoms][neighborList.len]
  Vector polarize;
  Vector relativeSolv;
  Vector solute;
} ForceField;

void readForceFieldFile(ForceField* forceField, char* forceFieldFile);
void forceFieldFree(ForceField* ff);
void vdwParameters(ForceField* ff, int nAtoms, int* atomClasses, Vector* neighborList);
void assignMultipoles(ForceField* forceField, REAL** multipoles, Vector* list12, Vector* list13, int* atomClasses, int nAtoms);

#endif //FORCEFIELDREADER_H