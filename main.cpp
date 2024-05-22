#include "Simulation.h"

int main(int argc, char* argv[])
{
    Simulation sim = Simulation();
    
    if(argc == 1)
    {
        // setup simulation parameters (temperature, nsteps, timestep, input files, output file, etc.)
        cout << "Reading input file" << endl;
        sim.readParameterFile("input.txt");
    }
    else
    {
        cout << "Reading " << argv[1] << endl;
        sim.readParameterFile(argv[1]);
    }
    
    // run simulation
    sim.run();
    
    return 0;
}