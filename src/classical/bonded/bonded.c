//
// Created by matthew-speranza on 8/20/24.
//

#include "../include/bonded.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/**
 * 3D vector subtraction.
 * @param a
 * @param b
 * @param c storage for the result
 */
void sub(REAL* a, REAL* b, REAL* c) {
  c[0] = a[0] - b[0];
  c[1] = a[1] - b[1];
  c[2] = a[2] - b[2];
}
void scale(REAL* a, REAL b, REAL* c) {
  c[0] = a[0] * b;
  c[1] = a[1] * b;
  c[2] = a[2] * b;
}
REAL length(REAL* a) {
  return sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
}
void normalize(REAL* a, REAL* b) {
  return scale(a, 1.0/length(a), b);
}
void add(REAL* a, REAL* b, REAL* c) {
  c[0] = a[0] + b[0];
  c[1] = a[1] + b[1];
  c[2] = a[2] + b[2];

}
REAL dot(REAL* a, REAL* b) {
  REAL dotVal = 0.0;
  for(int i = 0; i < 3; i++) {
    dotVal += a[i] * b[i];
  }
  return dotVal;
}

// This is identical to FFX version of these methods.
void rotateMultipole(System* system, int i) {
  // switch on multipole type
  int* frame = system->frameDef[i];
  enum MultipoleFrameDef frameDef = frame[4];
  REAL* localFrameMultipole = system->multipoles[i];
  REAL frameCoords[3][3];
  for(int j = 1; j < 4; j++) { // Ignore the atom itself
    if(frame[j] != -1) {
      int atomID = frame[j]*3;
      frameCoords[j-1][0] = system->X[atomID];
      frameCoords[j-1][1] = system->X[atomID+1];
      frameCoords[j-1][2] = system->X[atomID+2];
    }
  }
  REAL localOrigin[3] = {system->X[i*3], system->X[i*3+1], system->X[i*3+2]};
  REAL xAxis[3];
  REAL yAxis[3];
  REAL zAxis[3];

  // Use the identity matrix as the default rotation matrix
  REAL rotMat[3][3];
  rotMat[0][0] = 1.0;
  rotMat[1][0] = 0.0;
  rotMat[2][0] = 0.0;
  rotMat[0][1] = 0.0;
  rotMat[1][1] = 1.0;
  rotMat[2][1] = 0.0;
  rotMat[0][2] = 0.0;
  rotMat[1][2] = 0.0;
  rotMat[2][2] = 1.0;

  // Gen rotation matrix from frame def
  switch (frameDef) {
    case NONE:
      break;
    case BISECTOR:
      zAxis[0] = frameCoords[0][0];
      zAxis[1] = frameCoords[0][1];
      zAxis[2] = frameCoords[0][2];
      xAxis[0] = frameCoords[1][0];
      xAxis[1] = frameCoords[1][1];
      xAxis[2] = frameCoords[1][2];
      sub(zAxis, localOrigin, zAxis);
      normalize(zAxis, zAxis);
      sub(xAxis, localOrigin, xAxis);
      normalize(xAxis, xAxis);
      add(xAxis, zAxis, zAxis);
      normalize(zAxis, zAxis);
      rotMat[0][2] = zAxis[0];
      rotMat[1][2] = zAxis[1];
      rotMat[2][2] = zAxis[2];
      REAL dotVal = dot(xAxis, zAxis);
      scale(zAxis, dotVal, zAxis);
      sub(xAxis, zAxis, xAxis);
      normalize(xAxis, xAxis);
      break;
    case ZTHENBISECTOR:
      zAxis[0] = frameCoords[0][0];
      zAxis[1] = frameCoords[0][1];
      zAxis[2] = frameCoords[0][2];
      xAxis[0] = frameCoords[1][0];
      xAxis[1] = frameCoords[1][1];
      xAxis[2] = frameCoords[1][2];
      yAxis[0] = frameCoords[2][0];
      yAxis[1] = frameCoords[2][1];
      yAxis[2] = frameCoords[2][2];
      // Z-Axis
      sub(zAxis, localOrigin, zAxis);
      normalize(zAxis, zAxis);
      rotMat[0][2] = zAxis[0];
      rotMat[1][2] = zAxis[1];
      rotMat[2][2] = zAxis[2];
      sub(xAxis, localOrigin, xAxis);
      normalize(xAxis, xAxis);
      sub(yAxis, localOrigin, yAxis);
      normalize(yAxis, yAxis);
      // Sum the normalized vectors to the bisector atoms.
      add(xAxis, yAxis, xAxis);
      normalize(xAxis, xAxis);
      dotVal = dot(xAxis, zAxis);
      scale(zAxis, dotVal, zAxis);
      sub(xAxis, zAxis, xAxis);
      normalize(xAxis, xAxis);
      break;
    case ZONLY:
      zAxis[0] = frameCoords[0][0];
      zAxis[1] = frameCoords[0][1];
      zAxis[2] = frameCoords[0][2];
      sub(zAxis, localOrigin, zAxis);
      normalize(zAxis, zAxis);
      rotMat[0][2] = zAxis[0];
      rotMat[1][2] = zAxis[1];
      rotMat[2][2] = zAxis[2];
      // X-Axis: initially assume its along the global X-axis.
      xAxis[0] = 1.0;
      xAxis[1] = 0.0;
      xAxis[2] = 0.0;
      // If the Z-axis is close to the global X-axis,
      // assume an X-axis along the global Y-axis.
      dotVal = rotMat[0][2];
      if (fabs(dotVal) > 0.866) {
        xAxis[0] = 0.0;
        xAxis[1] = 1.0;
        dotVal = rotMat[1][2];
      }
      scale(zAxis, dotVal, zAxis);
      sub(xAxis, zAxis, xAxis);
      normalize(xAxis, xAxis);
      break;
    // 3-Fold frame rotation matrix elements for Z- and X-axes.
    case THREEFOLD:
      zAxis[0] = frameCoords[0][0];
      zAxis[1] = frameCoords[0][1];
      zAxis[2] = frameCoords[0][2];
      sub(zAxis, localOrigin, zAxis);
      normalize(zAxis, zAxis);
      xAxis[0] = frameCoords[1][0];
      xAxis[1] = frameCoords[1][1];
      xAxis[2] = frameCoords[1][2];
      sub(xAxis, localOrigin, xAxis);
      normalize(xAxis, xAxis);
      yAxis[0] = frameCoords[2][0];
      yAxis[1] = frameCoords[2][1];
      yAxis[2] = frameCoords[2][2];
      sub(yAxis, localOrigin, yAxis);
      normalize(yAxis, yAxis);
      add(zAxis, xAxis, zAxis);
      add(zAxis, yAxis, zAxis);
      normalize(zAxis, zAxis);
      rotMat[0][2] = zAxis[0];
      rotMat[1][2] = zAxis[1];
      rotMat[2][2] = zAxis[2];
      dotVal = dot(xAxis, zAxis);
      scale(zAxis, dotVal, zAxis);
      sub(xAxis, zAxis, xAxis);
      normalize(xAxis, xAxis);
      break;
    case ZTHENX:
      zAxis[0] = frameCoords[0][0];
      zAxis[1] = frameCoords[0][1];
      zAxis[2] = frameCoords[0][2];
      xAxis[0] = frameCoords[1][0];
      xAxis[1] = frameCoords[1][1];
      xAxis[2] = frameCoords[1][2];
      sub(zAxis, localOrigin, zAxis);
      normalize(zAxis, zAxis);
      rotMat[0][2] = zAxis[0];
      rotMat[1][2] = zAxis[1];
      rotMat[2][2] = zAxis[2];
      sub(xAxis, localOrigin, xAxis);
      dotVal = dot(xAxis, zAxis);
      scale(zAxis, dotVal, zAxis);
      sub(xAxis, zAxis, xAxis);
      normalize(xAxis, xAxis);
      break;
    default:
      printf("Couldn't recognize multipole frame definition type!\n");
  }
  // Set the X elements.
  rotMat[0][0] = xAxis[0];
  rotMat[1][0] = xAxis[1];
  rotMat[2][0] = xAxis[2];
  // Set the Y elements.
  rotMat[0][1] = rotMat[2][0] * rotMat[1][2] - rotMat[1][0] * rotMat[2][2];
  rotMat[1][1] = rotMat[0][0] * rotMat[2][2] - rotMat[2][0] * rotMat[0][2];
  rotMat[2][1] = rotMat[1][0] * rotMat[0][2] - rotMat[0][0] * rotMat[1][2];
  // Rotate the local frame multipole into the global frame.
  REAL r00 = rotMat[0][0];
  REAL r01 = rotMat[0][1];
  REAL r02 = rotMat[0][2];
  REAL r10 = rotMat[1][0];
  REAL r11 = rotMat[1][1];
  REAL r12 = rotMat[1][2];
  REAL r20 = rotMat[2][0];
  REAL r21 = rotMat[2][1];
  REAL r22 = rotMat[2][2];
  // Rotate the permanent dipole.
  REAL dx = localFrameMultipole[1];
  REAL dy = localFrameMultipole[2];
  REAL dz = localFrameMultipole[3];
  REAL* gMpole = system->rotatedMpoles[i];
  gMpole[0] = localFrameMultipole[0];
  gMpole[1] = r00 * dx + r01 * dy + r02 * dz;
  gMpole[2] = r10 * dx + r11 * dy + r12 * dz;
  gMpole[3] = r20 * dx + r21 * dy + r22 * dz;

  REAL qxx = localFrameMultipole[4] * 3;
  REAL qyy = localFrameMultipole[5] * 3;
  REAL qzz = localFrameMultipole[6] * 3;
  REAL qxy = localFrameMultipole[7] * 3;
  REAL qxz = localFrameMultipole[8] * 3;
  REAL qyz = localFrameMultipole[9] * 3;

  gMpole[4] = r00 * (r00 * qxx + r01 * qxy + r02 * qxz)
      + r01 * (r00 * qxy + r01 * qyy + r02 * qyz)
      + r02 * (r00 * qxz + r01 * qyz + r02 * qzz);
  gMpole[4] /= 3;

  gMpole[7] = r00 * (r10 * qxx + r11 * qxy + r12 * qxz)
      + r01 * (r10 * qxy + r11 * qyy + r12 * qyz)
      + r02 * (r10 * qxz + r11 * qyz + r12 * qzz);
  gMpole[7] /= 3;

  gMpole[8] = r00 * (r20 * qxx + r21 * qxy + r22 * qxz)
      + r01 * (r20 * qxy + r21 * qyy + r22 * qyz)
      + r02 * (r20 * qxz + r21 * qyz + r22 * qzz);
  gMpole[8] /= 3;

  gMpole[5] = r10 * (r10 * qxx + r11 * qxy + r12 * qxz)
      + r11 * (r10 * qxy + r11 * qyy + r12 * qyz)
      + r12 * (r10 * qxz + r11 * qyz + r12 * qzz);
  gMpole[5] /= 3;

  gMpole[9] = r10 * (r20 * qxx + r21 * qxy + r22 * qxz)
      + r11 * (r20 * qxy + r21 * qyy + r22 * qyz)
      + r12 * (r20 * qxz + r21 * qyz + r22 * qzz);
  gMpole[9] /= 3;

  gMpole[6] = r20 * (r20 * qxx + r21 * qxy + r22 * qxz)
      + r21 * (r20 * qxy + r21 * qyy + r22 * qyz)
      + r22 * (r20 * qxz + r21 * qyz + r22 * qzz);
  gMpole[6] /= 3;
};

void bondedTerms(System* system, int i) {
  //Compute and accumulate potentials and add grad into forces vector
  //bond(system, i)
  //piTorsion(system, i)
  //strBend(system, i)
  //angle(system, i)
  //oop(system, i)
  //torsion(system, i)
  //torsionTorsion(system, i)
}

void bondedLoop(System* system) {
  for(int i = 0; i < system->nAtoms; i++) {
    rotateMultipole(system, i);
    //bondedTerms(system, i);
  }
}
