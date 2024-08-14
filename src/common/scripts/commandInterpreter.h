#ifndef MOLECULAR_DYNAMICS_C_COMMANDINTERPRETER_H
#define MOLECULAR_DYNAMICS_C_COMMANDINTERPRETER_H

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
 void printLogo();
#endif //MOLECULAR_DYNAMICS_C_COMMANDINTERPRETER_H
