#ifndef MISCSTRUCTS_H
#define MISCSTRUCTS_H

#include "Coordinate.h"

struct Bond
{
	int i, j, bondTypeIndex;
    double creationTime;
    double destructionTime;
    
    Bond(int a, int b, int c) : i(a), j(b), bondTypeIndex(c)
    {
        creationTime = 0.0;
    }
    
    Bond(int a, int b, int c, double d) : i(a), j(b), bondTypeIndex(c), creationTime(d)
    {
    }
    
    Bond(int a, int b, int c, double d, double e) : i(a), j(b), bondTypeIndex(c), creationTime(d), destructionTime(e)
    {
    }
    
    bool operator<(const Bond& a) const
    {
        return (this->i < a.i);
    }
};

struct HalfBond
{
    int j, bondTag;
    
    HalfBond(int a, int b) : j{a}, bondTag{b}
    {
    }
};

struct Angle
{
	int i, j, k, angleTypeIndex;
    
    Angle(int a, int b, int c, int d) : i{a}, j{b}, k{c}, angleTypeIndex{d}
    {
    }
};

struct TemporaryBond
{
	int i, j;
	struct Coordinate force;
};

struct TemporaryBondAndDistance
{
    int i, j;
    struct Coordinate force;
    struct Coordinate distance;
};

struct TemporaryBondMag
{
    int i, j;
    double forceMag;
};

struct TemporaryForce
{
	int i;
	struct Coordinate force;
};

#endif