# Molecular Dynamics in C
No-one needs to see more production MD engines. This one is for learning.

This is the C version of this code base. I will try in the future to expand to other languages. Practice makes perfect.

## Usage
This project requires the clang compiler and cmake.

To build: 
- Create a build directory in the same directory as this readme file
- Execute "cmake -DCMAKE_BUILD_TYPE=Release ../" (change Release to Debug optionally)
- Then execute "make ."
- The executable "molecular_dynamics_c" should be in your build directory (move to bin or wherever you want)
- The executables "commonTest", "classicalTest", and "quantumTest" should also appear (for checking if tests pass)

To run static analysis:
- Execute "scan-build-{VERSION} cmake ../" from build directory

To Run:
- Have a valid structure
- Make a key file (see keyReader.h for keywords or grab one from examples/keyExamples/)
- Input four arguments "molecular_dynamics_c "command" "structure" "key" (or "molecular_dynamics_c help" for help menu)

## My Goals (Theory)
- Explain all the theory and reference heavily so this is educational for other people.
  
### Classical
- See what the simplest AMOEBA engine possible could look like.
- Do NPT, NVT, NVE simulations.
- Perform perturbative FE calculations via MBAR/BAR.
- Write a PB solver.
  
### Quantum
- Learn about foundational quantum chemistry algorithms (Hartree-Fock & Post) and theory following: https://github.com/jjgoings/McMurchie-Davidson
- Write the simplest possible DFT engine that can produce multipole components for AMOEBA. 

## My Goals (Programming)
- Learn the c programming language (no external libraries except perhaps OpenMP and CUDA)
- Write code that can be very simply parallelized on CPU and GPU by planning ahead.
