// Author(s): Matthew Speranza
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "../include/commandInterpreter.h"

#include <energy.h>

#include "../include/xyz.h"
#include "../include/keyReader.h"
#include "../include/neighborList.h"

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
    // Get structure file extension and read it in
    System* system = malloc(sizeof(System));
    if(system == NULL) {
        printf("malloc() failed to allocate memory in systemCreate()!");
        exit(1);
    }
    char* sExt = getFileExtension(structureFile, 3);
    assert(sExt != NULL);
    if(strcasecmp(sExt, supportedStructureExtensions[0]) == 0) { // xyz
        printf("Reading structure file: %s\n", structureFile);
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
        printf("Reading key file: %s\n", keyFile);
        readKeyFile(system, keyFile);
    } else {
        printf("Unsupported key file extension: %s\n", kExt);
        exit(1);
    }
    free(kExt);

    // Neighbors, Bonded terms/lists/params, atom classes
    buildLists(system);

    // Sets default settings if not already set & checks for complete system
    //setDefaults(system);
    // Multipoles to atoms and other stuff
    //finalizeForceField();
    // VdW Parameters (loops over neighbor lists created in last step)
    vdwParameters(system->forceField, system->nAtoms, system->atomClasses, system->verletList);

    return system;
}

/**
 * Frees all memory assiciated with a system.
 * @param system system to have all of its memory freed
 */
void systemDestroy(System* system) {
    for(int i = 0; i < system->nAtoms; i++) {
        //free(system->multipoles[i]);
        free(system->atomNames[i]);
        vectorBackingFree(&system->list12[i]);
        vectorBackingFree(&system->list13[i]);
        vectorBackingFree(&system->list14[i]);
        vectorBackingFree(&system->verletList[i]);
    }
    free(system->atomTypes);
    free(system->atomClasses);
    free(system->multipoles);
    free(system->atomNames);
    free(system->list12);
    free(system->list13);
    free(system->list14);
    // loop over bonded term vector of int* to free
    free(system->verletList);
    free(system->protons);
    free(system->valence);
    //for(int i = 0; i < system->pmeGridspace[0]; i++) {
    //    for(int j = 0; j < system->pmeGridspace[1]; j++) {
    //        free(system->pmeGrid[i][j]);
    //    }
    //    free(system->pmeGrid[i]);
    //}
    //free(system->pmeGrid);
    free(system->pmeGridspace);
    forceFieldFree(system->forceField);
    //free(system->pmeGridFlat);
    //free(system->DOF);
    //free(system->DOFFrc);
    free(system->X);
    free(system->M);
    free(system->V);
    free(system->A);
    free(system->F);
    free(system->lambdas);
    //free(system->thetas);
    //free(system->thetaM);
    //free(system->thetaV);
    //free(system->thetaA);
    //free(system->thetaF);
    //free(system->activeLambdas);
    free(system->remark);
    free(system->forceFieldFile);
    vectorBackingFree(&system->patchFiles);
    //free(system->keyFileName);
    //free(system->threadIDs);
    free(system);
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
