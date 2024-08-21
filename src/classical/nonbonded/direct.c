// Author(s): Matthew Speranza
#include "../include/direct.h"

#include <math.h>
#include <stdio.h>

void directLoop(System* system) {
  printf("Direct Loop\n");
  for(int i = 0; i < system->nAtoms; i++) {
    int i3 = i*3;
    Vector list = system->verletList[i];
    for(int j = 0; j < list.size; j++) {
      int j3 = ((int*)list.array)[j]*3;
      // Calculate R


      // Electrostatics
      // Potential Energy & Spatial Derivatives
/*
       Van der Waals
       Potential Energies
       CHARMM rules
       Lennard-Jones 12-6, 1-4 Rmin scaling, 1-2 and 1-3 off, GEOMETRIC EPS rule
       AMOEBA rules
       Buffered 14-7, 1-4 full scaling, 1-2 and 1-3 off, HHG EPS rule
       Buffered 14-7

       Uij = lambda.i * lambda.j * eps.ij * t1 * t2
       t1 = (1+del)^(n-m) / (alpha + (rho + del)^(n-m))
       t2 = [(1+gam) / (alpha + rho^m + gam)] - 2

       n = repulsivePower, m = dispersivePower
       n = 14, m = 7, (lennard-jones n=12, m=6)

       eps.ij = well depth from ff = (HHG) 4.0*(eps.i * eps.j) / ((sqrt(eps.i) + sqrt(eps.j))^2)
                                  || (GEOMETRIC) sqrt(eps.i)*sqrt(eps.j)

       del = .07, gam = .12 (lennard-jones del=0, gam=0)

       alpha = ffALPHA * ? (1-lami) * (1-lamj) ?
       if (soft) { ffALPHA = .25, lambda = lam^beta } else { alpha = 0, lambda = 1 }
       beta = 1.0

       rho = Rij / Rmin.ij
       Rmin.ij = (HHG) 2.0 * (ri^3 + rj^3) / (ri^2 + rj^2)  || (GEOMETRIC) 2.0 * sqrt(ri) * sqrt(rj)
       ri&j = vdw radii of atoms i and j

       Special reduction rules (mainly for hydrogens)
       - Place an atom closer to the another atom
       - ReductionFactor * (hydrogen-x - rx) + rx
*/
      REAL dx = system->X[i3] - system->X[j3];
      REAL dy = system->X[i3+1] - system->X[j3+1];
      REAL dz = system->X[i3+2] - system->X[j3+1];
      REAL rij = sqrt(dx*dx + dy*dy + dz*dz);
      REAL lambdaI = system->lambdas[i];
      REAL lambdaJ = system->lambdas[j];
      REAL alpha = 0.0;
      REAL epsIJ = system->forceField->wellDepths[i][j];
      REAL rMin = system->forceField->rMin[i][j];
      REAL rho = rij / rMin;
      REAL n = system->forceField->vdwN;
      REAL m = system->forceField->vdwM;
      REAL delta = system->forceField->vdwDelta;
      REAL gamma = system->forceField->vdwGamma;
      REAL t1 = pow(1+delta, n-m) / (alpha + pow(rho + delta, n-m));
      REAL t2 = (1+gamma) / (alpha + pow(rho, m) + gamma) - 2;
      REAL uij = lambdaI * lambdaJ * epsIJ * t1 * t2;
      //printf("U.vdw %d %d: %lf\n", i, j, uij);
      // Spatial Derivatives
      // Lambda Derivatives
    }
  }
};
