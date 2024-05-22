#include "SimulationBox.h"

SimulationBox::SimulationBox()
{
    boxLengths = Coordinate {1.0, 1.0, 1.0};
    periodicX = true;
    periodicY = true;
    periodicZ = true;
}

SimulationBox::SimulationBox(Coordinate b, bool px, bool py, bool pz)
{
    boxLengths = b;
    
    periodicX = px;
    periodicY = py;
    periodicZ = pz;
}

double SimulationBox::getBoxLength(int index)
{
    if(index == 0)
        return boxLengths.x;
    else if(index == 1)
        return boxLengths.y;
    else if(index == 2)
        return boxLengths.z;
    else
        throw std::runtime_error("SimulationBox does not have an index of " + std::to_string(index));
    
    return 0.0;
}

void SimulationBox::setBoxLength(int index, double newVal)
{
    if(index == 0)
        boxLengths.x = newVal;
    else if(index == 1)
        boxLengths.y = newVal;
    else if(index == 2)
        boxLengths.z = newVal;
}

double SimulationBox::periodicWrap(double pos, double boxl)
{
    return pos - boxl * round(pos/boxl);
}

Coordinate SimulationBox::periodicWrap(Coordinate a)
{
    if(periodicX)
        a.x = periodicWrap(a.x, boxLengths.x);
    if(periodicY)
        a.y = periodicWrap(a.y, boxLengths.y);
    if(periodicZ)
        a.z = periodicWrap(a.z, boxLengths.z);
    
    return a;
}

struct Coordinate SimulationBox::calcDisplacement(Coordinate a, Coordinate b)
{
    Coordinate c = a - b;
    
    c = periodicWrap(c);
    
    return c;
}

double SimulationBox::calcDistance(Coordinate a, Coordinate b)
{
    Coordinate c = calcDisplacement(a, b);
    
    return c.getMagnitude();
}