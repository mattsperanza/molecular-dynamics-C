// Author(s): Matthew Speranza
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>

#include "../include/commandInterpreter.h"

#include <energy.h>
#include <omp.h>

#include "../include/xyz.h"
#include "../include/keyReader.h"
#include "../include/neighborList.h"
#include "timer.h"

int nSupStructExt = 3;
char* supportedStructureExtensions[3] = {"xyz", "arc", "pdb"};
int nSupKeyExt = 2;
char* supportedKeyExtensions[2] = {"key", "properties"};

/**
 * @param argc number of arguments, must be four
 * @param argv array of arguments - expects [${PATH}/molecular_dynamics_c, command, [supported structure file], key file]
 */
void commandInterpreter(int argc, char *argv[]){
    printLogo();
    printf("Command line input: \"");
    // Ignore path up do md
    char *firstArg = argv[0];
    int lastSeparator = 0;
    for(int i = 0; i < strlen(firstArg); i++){
        if(firstArg[i] == '/' || firstArg[i] == '\\'){
            lastSeparator = i;
        }
    }
    // Print each argument
    for(int i = lastSeparator + 1; i < strlen(firstArg); i++){
        printf("%c", firstArg[i]);
    }
    for (int i = 1; i < argc; i++) {
        printf(" %s", argv[i]);
    }
    printf("\"\n");

    if (argc < 2) {
        printf("No command provided. Exiting.\n");
        return;
    }
    char* command = argv[1];
    if(strcasecmp(command, "help") == 0) {
        printf("Entering help menu.");
    } else if(strcasecmp(command, "energy") == 0 && argc == 4) {
        printf("Preparing to calculate the energy of the system.\n");
        System* system = systemCreate(argv[2], argv[3]);
        energy(system);
        systemDestroy(system);
    } else if(strcasecmp(command, "dynamics") == 0 && argc == 4) {
        printf("Preparing to run molecular dynamics on the system.\n");
        System* system = systemCreate(argv[2], argv[3]);
        //dynamics(system); // calls energy many times
        systemDestroy(system);
    } else if (argc != 4){
        printf("Program expects 3 arguments in addition to command if help isn't requested!\n");
        printf("Required format: \"[$COMMAND_PATH, supported command, supported structure file, key file");
        printf("Supported commands: ");
        printSupportedCommands();
        printf("Supported structure files: ");
        printSupportedStructureFiles();
    }
}

char* getFileExtension(char* fileName, int extForceLen) {
    int lastDotIndex = 0;
    for(int i = strlen(fileName)-1; i >= 0; i--) {
       char letter = fileName[i];
       if(letter == '.') {
           lastDotIndex = i;
           break;
       }
       if(i == 0) {
           printf("Failed to read file extension from: %s\n", fileName);
           exit(1);
       }
    }
    // add one for end char
    int extLen = extForceLen != -1 ? extForceLen + 1: strlen(fileName) - (lastDotIndex+1) + 1;
    if(extLen > strlen(fileName) - (lastDotIndex+1) + 1) {
        printf("Invalid extension length given in getFileExtension!");
        exit(1);
    }
    char* ext = malloc(sizeof(char)*extLen+1);
    if(ext == NULL) {
        printf("malloc() failed to allocate memory in getFileExtension()!");
        exit(1);
    }
    for(int i = 0; i < extLen; i++) {
        ext[i] = fileName[lastDotIndex+1+i];
    }
    return ext;
}

System* systemCreate(char* structureFile, char* keyFile) {
    double sysStart, sysStop;
    // Get structure file extension and read it in
    sysStart = omp_get_wtime();
    System* system = malloc(sizeof(System));
    if(system == NULL) {
        printf("malloc() failed to allocate memory in systemCreate()!");
        exit(1);
    }
    char* sExt = getFileExtension(structureFile, 3);
    assert(sExt != NULL);
    if(strcasecmp(sExt, supportedStructureExtensions[0]) == 0) { // xyz
        readXYZ(system, structureFile);
    } else if (strcasecmp(sExt, supportedStructureExtensions[1]) == 0){ // arc -> extended xyz
        printf("Reading structure file: %s\n", structureFile);
        //readARC(system, structureFile);
    } else if (strcasecmp(sExt, supportedStructureExtensions[2]) == 0) { // pdb
        printf("PDB reader has not been implemented yet!\n");
        exit(1);
    } else {
        printf("Unsupported structure file extension: %s\n", sExt);
        printf("Supported extensions: ");
        printSupportedStructureFiles();
        exit(1);
    }
    free(sExt);

    // Key file reader - also reads force field file
    char* kExt = getFileExtension(keyFile,-1);
    assert(kExt != NULL);
    if(strcasecmp(kExt, supportedKeyExtensions[0]) == 0 || strcasecmp(kExt, supportedKeyExtensions[1]) == 0) { // key file
        readKeyFile(system, keyFile);
    } else {
        printf("Unsupported key file extension: %s\n", kExt);
        exit(1);
    }
    free(kExt);
    sysStop = omp_get_wtime();
    printClock(sysStart, sysStop, "read system files");

    // Sets default settings if not already set & checks for complete system
    setDefaults(system);
    checkSystem(system);

    // Neighbors, Bonded terms/lists/params, atom classes
    double listStart, listStop;
    listStart = omp_get_wtime();
    buildLists(system);
    listStop = omp_get_wtime();
    printClock(listStart, listStop, "assign bonded lists and params");

    // Matches atoms with multipole def in FF file and defines frame atom indices with sign of atom type in ff file
    // Allocs mem for multipoles
    double assignStart, assignStop;
    assignStart = omp_get_wtime();
    assignMultipoles(system->forceField, &system->multipoles, &system->rotatedMpoles,
        &system->frameDef, system->list12, system->list13, system->atomTypes, system->nAtoms);
    // VdW Parameters (loops over neighbor lists created in neighbor list step) - needs to get called every neighbor list rebuild
    vdwParameters(system->forceField, system->nAtoms, system->atomClasses, system->verletList);
    assignStop = omp_get_wtime();
    printClock(assignStart, assignStop, "assign non-bonded params");

    return system;
}

void checkSystem(const System* system) {
    if(system->X == NULL) {
        printf("Coordinates have not been read in yet!\n");
        exit(1);
    }
    if(system->M == NULL) {
        printf("Masses are still NULL!\n");
        exit(1);
    }
    if(system->V == NULL) {
        printf("Velocities are still NULL!\n");
        exit(1);
    }
    if(system->A == NULL) {
        printf("Accelerations are still NULL!\n");
        exit(1);
    }
    if(system->F == NULL) {
        printf("Forces are still NULL!\n");
        exit(1);
    }
    if(system->lambdas == NULL) {
        printf("Lambdas are still NULL!\n");
        exit(1);
    }
    for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
            if(system->boxDim[i][j] == -1.0) {
                printf("One of the box dimensions is still -1.0!\n");
                exit(1);
            }
            if(i == j && system->boxDim[i][j] == 0.0){
                printf("Box is not 3D!");
                exit(1);
            }
        }
    }
}

void setDefaults(System* system) {
    printf("Setting system defaults:\n");
    if(system->temperature == 0 || system->temperature < 0) {
        printf(" - Defaulting to 298.15 Kelvin for temperature.\n");
        system->temperature = 298.15;
    }
    if(system->dtFemto == 0 || system->dtFemto < 0) {
        printf(" - Defaulting to 1 femtoseconds per timestep.\n");
        system->dtFemto = 1.0;
        system->dtAtto = 1000;
    }
    if(system->printThermoEvery == 0 || system->printThermoEvery < 0) {
        printf(" - Defaulting to printing energy every 1000 femtoseconds.\n");
        system->printThermoEvery = 1000;
    }
    if(system->printRestartEvery == 0 || system->printRestartEvery < 0) {
        printf(" - Defaulting to printing restart file every 1000 femtoseconds.\n");
        system->printRestartEvery = 1000;
    }
    if(system->ewaldAlpha == 0.0) {
        printf(" - EWALD SUMMATION IS TURNED OFF!\n");
        printf(" - Note that ewald alpha is 0.0 by default and must be set (to ~0.545) for Ewald summation typically.\n");
    }
    if(system->ewaldOrder == 0 || system->ewaldOrder < 0) {
        printf(" - Defaulting to 5 for Ewald order.\n");
        system->ewaldOrder = 5;
    }
    if(system->pmeGridspace[0] == 0 || system->pmeGridspace[0] < 0) {
        printf(" - Defaulting to 32 gridpoints in x for PME.\n");
        system->pmeGridspace[0] = 32;
    }
    if(system->pmeGridspace[1] == 0 || system->pmeGridspace[1] < 0) {
        printf(" - Defaulting to 32 gridpoints in y for PME.\n");
        system->pmeGridspace[1] = 32;
    }
    if(system->pmeGridspace[2] == 0 || system->pmeGridspace[2] < 0) {
        printf(" - Defaulting to 32 gridpoints in z for PME.\n");
        system->pmeGridspace[2] = 32;
    }
    if(system->ewaldCutoff == 0.0 || system->ewaldCutoff < 0) {
        printf(" - Defaulting to 7.0 Angstroms for Ewald cutoff.\n");
        system->ewaldCutoff = 7.0;
    }
    if(system->vdwCutoff == 0.0 || system->vdwCutoff < 0) {
        printf(" - Defaulting to 12.0 Angstroms for VdW cutoff.\n");
        system->vdwCutoff = 12.0;
    }
    if(system->vdwTaper == 0.0 || system->vdwTaper < 0) {
        printf(" - Defaulting to 0.9 for VdW tapering.\n");
        system->vdwTaper = 0.9;
    }
    if(system->realspaceBuffer == 0.0 || system->realspaceBuffer < 0) {
        printf(" - Defaulting to 2 Angstroms for real space buffer.\n");
        system->realspaceBuffer = 2.0;
    }
}

/**
 * Frees all memory assiciated with a system.
 * @param system system to have all of its memory freed
 */
void systemDestroy(System* system) {
    // This list of frees is up to date-ish but couldn't bother to really think about this since it's at the end anyway
    for(int i = 0; i < system->nAtoms; i++) {
        //free(system->multipoles[i]);
        //free(system->atomNames[i]);
        //vectorBackingFree(&system->list12[i]);
        //vectorBackingFree(&system->list13[i]);
        //vectorBackingFree(&system->list14[i]);
        //vectorBackingFree(&system->verletList[i]);
    }
    //free(system->atomTypes);
    //free(system->atomClasses);
    //free(system->multipoles);
    //free(system->atomNames);
    //free(system->list12);
    //free(system->list13);
    //free(system->list14);
    // loop over bonded term vector of int* to free
    //free(system->verletList);
    //free(system->protons);
    //free(system->valence);
    //for(int i = 0; i < system->pmeGridspace[0]; i++) {
    //    for(int j = 0; j < system->pmeGridspace[1]; j++) {
    //        free(system->pmeGrid[i][j]);
    //    }
    //    free(system->pmeGrid[i]);
    //}
    //free(system->pmeGrid);
    //free(system->pmeGridspace);
    //forceFieldFree(system->forceField);
    //free(system->pmeGridFlat);
    //free(system->DOF);
    //free(system->DOFFrc);
    //free(system->X);
    //free(system->M);
    //free(system->V);
    //free(system->A);
    //free(system->F);
    //free(system->lambdas);
    //free(system->thetas);
    //free(system->thetaM);
    //free(system->thetaV);
    //free(system->thetaA);
    //free(system->thetaF);
    //free(system->activeLambdas);
    //free(system->remark);
    //free(system->forceFieldFile);
    //vectorBackingFree(&system->patchFiles);
    //free(system->keyFileName);
    //free(system->threadIDs);
    //free(system);
}

void printSupportedCommands() {
    int nCommands = 3;
    char* commands[3] = {"help", "energy", "dynamics"};
    printf("[ ");
    for(int i = 0; i < nCommands; i++) {
        printf("%s ", commands[i]);
    }
    printf("]");
}

void printSupportedStructureFiles() {
    int nStructureTypes = 3;
    char* structures[3] = {"*.pdb", "*.xyz", "*.arc"};
    printf("[ ");
    for(int i = 0; i < nStructureTypes; i++) {
        printf("%s ", structures[i]);
    }
    printf("]");
}

void printLogo() {
    char *logoString;
    logoString = "------------------------------------------------------------------------------------------\n"
                 "|   ███╗   ███╗ ██████╗ ██╗     ███████╗ ██████╗██╗   ██╗██╗      █████╗ ██████╗         |\n"
                 "|   ████╗ ████║██╔═══██╗██║     ██╔════╝██╔════╝██║   ██║██║     ██╔══██╗██╔══██╗        |\n"
                 "|   ██╔████╔██║██║   ██║██║     █████╗  ██║     ██║   ██║██║     ███████║██████╔╝        |\n"
                 "|   ██║╚██╔╝██║██║   ██║██║     ██╔══╝  ██║     ██║   ██║██║     ██╔══██║██╔══██╗        |\n"
                 "|   ██║ ╚═╝ ██║╚██████╔╝███████╗███████╗╚██████╗╚██████╔╝███████╗██║  ██║██║  ██║        |\n"
                 "|   ╚═╝     ╚═╝ ╚═════╝ ╚══════╝╚══════╝ ╚═════╝ ╚═════╝ ╚══════╝╚═╝  ╚═╝╚═╝  ╚═╝        |\n"
                 "|                    ██████╗ ██╗   ██╗███╗   ██╗ █████╗ ███╗   ███╗██╗ ██████╗███████╗   |\n"
                 "|                    ██╔══██╗╚██╗ ██╔╝████╗  ██║██╔══██╗████╗ ████║██║██╔════╝██╔════╝   |\n"
                 "|       ███████╗     ██║  ██║ ╚████╔╝ ██╔██╗ ██║███████║██╔████╔██║██║██║     ███████╗   |\n"
                 "|       ╚══════╝     ██║  ██║  ╚██╔╝  ██║╚██╗██║██╔══██║██║╚██╔╝██║██║██║     ╚════██║   |\n"
                 "|                    ██████╔╝   ██║   ██║ ╚████║██║  ██║██║ ╚═╝ ██║██║╚██████╗███████║   |\n"
                 "|                    ╚═════╝    ╚═╝   ╚═╝  ╚═══╝╚═╝  ╚═╝╚═╝     ╚═╝╚═╝ ╚═════╝╚══════╝   |\n"
                 "|                                                                                        |\n"
                 "|                                                                                        |\n"
                 "|             ██████╗      ███████╗██████╗ ██╗████████╗██╗ ██████╗ ███╗   ██╗            |\n"
                 "|            ██╔════╝      ██╔════╝██╔══██╗██║╚══██╔══╝██║██╔═══██╗████╗  ██║            |\n"
                 "|            ██║     █████╗█████╗  ██║  ██║██║   ██║   ██║██║   ██║██╔██╗ ██║            |\n"
                 "|            ██║     ╚════╝██╔══╝  ██║  ██║██║   ██║   ██║██║   ██║██║╚██╗██║            |\n"
                 "|            ╚██████╗      ███████╗██████╔╝██║   ██║   ██║╚██████╔╝██║ ╚████║            |\n"
                 "|             ╚═════╝      ╚══════╝╚═════╝ ╚═╝   ╚═╝   ╚═╝ ╚═════╝ ╚═╝  ╚═══╝            |\n"
                 "|   Written by Matt.                                                                     |\n"
                 "------------------------------------------------------------------------------------------\n";
    printf("%s\n\n", logoString);
}
