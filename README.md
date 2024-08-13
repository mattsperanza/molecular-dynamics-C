# Molecular Dynamics in C
No-one needs to see more production MD engines. This one is for learning.

I want to write this in c first, then haskell, then zig to get really good understanding of numerics/struggles/theory/details of writing MD code.

Practice makes perfect.

## My Goals (Theory)
### Classical
- I want to see what the simplest AMOEBA engine possible could look like.
- I want to explain all of the theory and reference heavily so this is educational for other people.
- I want to do NPT, NVT, NVE simulations.
- I want to erform perturbative FE calculations via MBAR/BAR.
  
### Quantum
- I want to learn about foundational quantum chemistry algorithms (Hartree-Fock & Post) and theory following: https://github.com/jjgoings/McMurchie-Davidson
- I want to write the simplest possible DFT engine that can produce multipole components for AMOEBA. 

## My Goals (Programming)
- I want to learn the c programming language.
- I want to write code that can be very simply parallelized on CPU and GPU by planning ahead.
