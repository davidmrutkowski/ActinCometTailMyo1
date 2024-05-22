#include "SpringBondInfo.h"

using namespace std;

SpringBondInfo::SpringBondInfo()
{
	k = 0.0;
	eqDist = 0.0;
}

SpringBondInfo::SpringBondInfo(double newK, double newEqDist)
{
	k = newK;
	eqDist = newEqDist;
}

double SpringBondInfo::calcForce(double dist)
{
    //positive values are compression
	return -k * (dist - eqDist);
}

double SpringBondInfo::getEqDist()
{
    return eqDist;
}