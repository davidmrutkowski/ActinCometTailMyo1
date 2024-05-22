#ifndef SIMULATIONBOX_H
#define SIMULATIONBOX_H

#include <iostream>
#include <stdexcept>
#include "MiscStructs.h"

using namespace std;

class SimulationBox
{
	private:
		Coordinate boxLengths;
        bool periodicX, periodicY, periodicZ;
        
	public:
        SimulationBox();
		SimulationBox(Coordinate b, bool px, bool py, bool pz);
		
        double periodicWrap(double pos, double boxL);
        Coordinate periodicWrap(Coordinate a);
        struct Coordinate calcDisplacement(Coordinate c, Coordinate d);
        double calcDistance(Coordinate c, Coordinate d);
        double getBoxLength(int index);
        void setBoxLength(int index, double newVal);
};

#endif