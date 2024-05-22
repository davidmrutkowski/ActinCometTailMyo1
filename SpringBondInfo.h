#ifndef SPRINGBONDINFO_H
#define SPRINGBONDINFO_H

#include "BondInfo.h"

class SpringBondInfo:public BondInfo
{
	private:
		double k;
		double eqDist;
	public:
		SpringBondInfo();
		SpringBondInfo(double newK, double newEqDist);
		double calcForce(double dist);
        double getEqDist();
};

#endif