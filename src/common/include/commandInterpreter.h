// Author(s): Matthew Speranza
#ifndef MOLECULAR_DYNAMICS_C_COMMANDINTERPRETER_H
#define MOLECULAR_DYNAMICS_C_COMMANDINTERPRETER_H

#include "../system/system.h"

/**
 * Reads commands from user, parses supported files given, and moves flow to individual command handlers.
 * <hr>
 * Currently supported commands:
 * \Energy calculates the energy of the system
 * <p>
 * \Dynamics runs molecular dynamics on the system
 * <p>
 */
 void commandInterpreter(int argc, char *argv[]);
 void printSupportedCommands();
 void printSupportedStructureFiles();
 System* systemCreate(char* structureFileName, char* keyFileName);
 void systemDestroy(System* system); // I wanna move this to system.h but got linker errors
 char* getFileExtension(char* fileName, int extForceLen);
 void printLogo();
#endif //MOLECULAR_DYNAMICS_C_COMMANDINTERPRETER_H
