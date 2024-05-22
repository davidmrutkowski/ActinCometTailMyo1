#ifndef FILAMENT_H
#define FILAMENT_H

#include <vector>
#include <queue>
#include <iostream>
#include "Coordinate.h"
#include "ParticleInfo.h"

class Filament
{
	private:
		std::deque <int> p_tags;
		bool isDaughter;
		int motherTag;
        bool cappedBarbedEnd;
        bool cappedPointedEnd;
	public:
		Filament();
		
		int getNumParticles();
		int getTagAtIndex(int pos);
		
		int addParticleFront(int newTag);
		int addParticleBack(int newTag);
		void removeParticleFront();
		void removeParticleBack();
		
		struct Coordinate calcBendingPotential(double k, ParticleInfo &pinfo, struct TemporaryBondMag *maxBond, bool bendingInStress);

		struct Coordinate getFrontDirection(ParticleInfo pinfo);
		
		static struct Coordinate matrixProduct(double basis[3][3], struct Coordinate a);
		
        void setMotherTag(int newMotherTag);
		int getMotherTag();
        
		void setIsDaughter(bool state);
		bool getIsDaughter();
		
        Filament breakFilamentAtPos(int pos);
        Filament removeBeadAtPos(int pos);
        Filament removeBeadByTag(int tag);
        
        int getPointedTag();
        int getBarbedTag();
        
        void setCappedBarbedEnd(bool newState);
        bool getCappedBarbedEnd();
        void setCappedPointedEnd(bool newState);
        bool getCappedPointedEnd();
        
        int getIndexOfTag(int tempTag);
        
        std::deque <int> getAllTags();
};

#endif