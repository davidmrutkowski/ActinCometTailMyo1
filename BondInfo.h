#ifndef BONDINFO_H
#define BONDINFO_H

class BondInfo
{
	public:
		virtual double calcForce(double dist)
		{ return 0.0; }
};

#endif