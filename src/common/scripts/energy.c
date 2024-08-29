// Author(s): Matthew Speranza

#include "../include/energy.h"
#include "../../classical/include/direct.h"
#include "../../classical/include/bonded.h"

void energy(System* system) {
  bondedLoop(system); // Rotates multipoles
  directLoop(system);
  //reciprocalLoop(system);
  //reciprocalPME(system);
  //scf(system);
};
