#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>
#include <map>
#include <queue>
#include <array>
#include <math.h>
#include "Coordinate.h"

class Particle
{
	private:
		std::array <double, 3> r;
		std::array <double, 3> force;
		double mass;
		int bin;
        int type;
		
	public:
		Particle();
		Particle(double rx, double ry, double rz);
		Particle(const Particle &p2);
		void setPosition(int index, double value);
		double getPosition(int index);
		void setForce(int index, double value);
		double getForce(int index);
		void zeroAllForces();
		void addForce(int index, double value);
		double getMass();
		void setBin(int newBin);
		int getBin();
        void setType(int newType);
        int getType();
};

#endif // USER_H