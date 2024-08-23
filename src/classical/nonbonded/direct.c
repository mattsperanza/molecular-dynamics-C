// Author(s): Matthew Speranza
#include "../include/direct.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void vdwTaperConstants(REAL a, REAL b, REAL c[6]) { // From MultiplicativeSwitch.java in FFX
  // a - taper from | b - taper to
  REAL a2 = a * a;
  REAL b2 = b * b;
  REAL denom = pow(b-a, 5.0);
  c[0] = b*b2*(b2 - 5.0*a*b + 10.0*a2) / denom;
  c[1] = -30.0*a2*b2/denom;
  c[2] = 30.0*b*a*(b+a)/denom;
  c[3] = -10.0 * (a2 + 4.0 * a * b + b2) / denom;
  c[4] = 15.0 * (a + b) / denom;
  c[5] = -6.0 / denom;
}

REAL vdwTaperFunction(REAL r, REAL c[6]) {
  REAL r2 = r * r;
  REAL r3 = r2 * r;
  REAL r4 = r2 * r2;
  REAL r5 = r4 * r;
  return c[5]*r5 + c[4]*r4 + c[3]*r3 + c[2]*r2 + c[1]*r + c[0];
}

REAL vdwTaperFunctionDerivative(REAL r, REAL c[6]) {
  REAL r2 = r * r;
  REAL r3 = r2 * r;
  REAL r4 = r2 * r2;
  return 5.0*c[5]*r4 + 4.0*c[4]*r3 + 3.0*c[3]*r2 + 2.0*c[2]*r + c[1];
}


void directLoop(System* system) {
  printf("Direct Loop\n");

  // VdW Params
  REAL taperConstants[6];
  REAL taperStart = system->realspaceCutoff * .9; // Default in FFX is .9
  vdwTaperConstants(taperStart, system->realspaceCutoff, taperConstants);
  REAL vdwE = 0.0;
  long vdwCount = 0;

  for(int i = 0; i < system->nAtoms; i++) {
    int i3 = i*3;
    int* list = system->verletList[i].array;
    REAL xiV = system->X[i3]; // Coord to use for vdw
    REAL yiV = system->X[i3+1];
    REAL ziV = system->X[i3+2];
    REAL redi = system->forceField->reductionFactors[i];
    if(redi != 0) { // Shift atom close to bonded heavy atom
      int heavy = ((int*)system->list12[i].array)[0]; // one heavy atom
      REAL heavyX = system->X[heavy*3];
      REAL heavyY = system->X[heavy*3+1];
      REAL heavyZ = system->X[heavy*3+2];
      xiV = redi * (xiV - heavyX) + heavyX;
      yiV = redi * (yiV - heavyY) + heavyY;
      ziV = redi * (ziV - heavyZ) + heavyZ;
    }
    REAL* vdwMask = calloc(sizeof(REAL), system->nAtoms);
    for(int j = 0; j < system->nAtoms; j++) {
      vdwMask[j] = 1.0;
    }
    // 1-2, 1-3 set to 0
    for(int jj = 0; jj < system->list12[i].size; jj++) {
      int j = ((int*)system->list12[i].array)[jj];
      vdwMask[j] = 0.0;
    }
    for(int jj = 0; jj < system->list13[i].size; jj++) {
      int j = ((int*)system->list13[i].array)[jj];
      vdwMask[j] = 0.0;
    }

    for(int jj = 0; jj < system->verletList[i].size; jj++) {
      int j = list[jj]; // Access system from this variable
      int j3 = j*3;
      // Electrostatics - see explanation in method
      // multipoleInteraction(system, i, j);


      // VDW References:
      //
      // ForceFieldX: https://github.com/SchniedersLab/forcefieldx
      //
      // Schnieders et al. The Structure, Thermodynamics, and Solubility of Organic Crystals
      // from Simulation with a Polarizable Force Field
      // https://pubs.acs.org/doi/epdf/10.1021/ct300035u
      //
      // CHARMM ALF soft-core: https://pubs.acs.org/doi/full/10.1021/acs.jpcb.6b09656
      //
      /* Van der Waals - Potential Energies
       * CHARMM rules:
       * Lennard-Jones 12-6, 1-4 Rmin scaling, 1-2 and 1-3 off, GEOMETRIC EPS rule
       * AMOEBA rules:
       * Buffered 14-7, 1-4 full scaling, 1-2 and 1-3 off, HHG EPS rule
       *
       * Uij = lambda.i * lambda.j * eps.ij * t1 * t2
       * t1 = (1+del)^(n-m) / ((rho + del)^(n-m))
       * t2 = [(1+gam) / (rho^m + gam)] - 2
       *
       * n = repulsivePower, m = dispersivePower
       * n = 14, m = 7, (lennard-jones n=12, m=6)
       *
       * eps.ij = well depth from ff
       * eps.ij (HHG) = 4.0*(eps.i * eps.j) / ((sqrt(eps.i) + sqrt(eps.j))^2
       * eps.ij (GEOMETRIC) = sqrt(eps.i)*sqrt(eps.j)
       *
       * del = .07, gam = .12 (lennard-jones del=0, gam=0)
       *
       * rho = Rij / Rmin.ij
       * Rmin.ij = (HHG) 2.0 * (ri^3 + rj^3) / (ri^2 + rj^2)  || (GEOMETRIC) 2.0 * sqrt(ri) * sqrt(rj)
       * ri&j = vdw radii of atoms i and j
       *
       * AMOEBA reduction rules (mainly for hydrogens)
       * - Place an atom closer to the another atom
       * - newHydrogenX = reductionFactor * (hydrogenX - heavyX) + heavyX
       */
      REAL xjV = system->X[j3];
      REAL yjV = system->X[j3+1];
      REAL zjV = system->X[j3+2];
      REAL redj = system->forceField->reductionFactors[j];
      if(redj != 0.0) {
        int heavy = ((int*)system->list12[j].array)[0];
        REAL heavyX = system->X[heavy*3];
        REAL heavyY = system->X[heavy*3+1];
        REAL heavyZ = system->X[heavy*3+2];
        xjV = redj * (xjV - heavyX) + heavyX;
        yjV = redj * (yjV - heavyY) + heavyY;
        zjV = redj * (zjV - heavyZ) + heavyZ;
      }
      REAL dx = pbc(xiV - xjV, system->boxDim[0][0]);
      REAL dy = pbc(yiV - yjV, system->boxDim[1][1]);
      REAL dz = pbc(ziV - zjV, system->boxDim[2][2]);
      REAL rijVdW = sqrt(dx*dx + dy*dy + dz*dz);
      REAL rMin = system->forceField->rMin[i][jj];
      REAL mask = vdwMask[j];
      if(mask <= 0.0 || rMin <= 0.0 || rijVdW > system->realspaceCutoff) {
        continue;
      }
      REAL lambdaI = system->lambdas[i];
      REAL lambdaJ = system->lambdas[j];
      REAL epsIJ = system->forceField->wellDepths[i][jj]; // these are relative to neighborlist
      REAL rho = rijVdW / rMin;
      REAL n = system->forceField->vdwN;
      REAL m = system->forceField->vdwM;
      REAL delta = system->forceField->vdwDelta;
      REAL gamma = system->forceField->vdwGamma;
      REAL t1 = pow(1+delta, n-m) / pow(rho + delta, n-m);
      REAL t2 = (1+gamma) / (pow(rho, m) + gamma) - 2;
      REAL taper = 1.0;
      if(rijVdW > system->realspaceCutoff - 1.2) {
        taper = vdwTaperFunction(rijVdW, taperConstants);
      }
      REAL uij = taper * mask * lambdaI * lambdaJ * epsIJ * t1 * t2;
      vdwE += uij;
      vdwCount++;

      /* VdW - Spatial Derivatives
       *
       *  Uij = lambda.i * lambda.j * eps.ij * t1 * t2
       *  t1 = (1+del)^(n-m) / ((rho + del)^(n-m))
       *  t2 = [(1+gam) / (rho^m + gam)] - 2
       *
       *
       *
       *
       *
       */

      /* VdW - Lambda Derivatives
       *
       *
       *
       *
       *
       */
    }
    free(vdwMask);
  }
  printf("VdW Count: %ld\n", vdwCount);
  printf("VdW Energy: %lf\n", vdwE);
};
