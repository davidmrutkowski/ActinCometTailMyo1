# ActinCometTailMyo1
Brownian dynamics code that simulates actin network growth (actin branching and filament elongation) around a nucleator bead under the effect of myosin-I. This code was developed for the simulations in Xu, Rutkowski, Rebowski, Boczkowska, Pollard, Dominguez, Vavylonis, and Ostap, bioRiv, https://doi.org/10.1101/2024.02.09.579714

Compilation using g++: g++ -O3 -fopenmp -o main main.cpp *.cpp

Input parameter file (input.txt): This file lists various parameters used by the simulation.
Command to run code: ./main input.txt

Restart from previously run simulation: Add the following lines to input.txt: 
positions last###.xyz 
bonds last###.bnd
angles last###.xyz
Change the outputName parameter to a new name to keep from overwriting previous output.
