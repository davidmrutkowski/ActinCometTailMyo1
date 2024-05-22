#ifndef GRID_H
#define GRID_H

#include <vector>
#include <tuple>
#include <algorithm>
#include <iostream>
#include "Coordinate.h"

class Grid
{
	private:
		std::vector<std::vector<int>> grid;
        std::vector<int> tags;
		int numBinX, numBinY, numBinZ;
		Coordinate gridSize;
		Coordinate binSize;
		bool periodicX, periodicY, periodicZ;
	public:
        Grid();
        Grid(Coordinate b, bool px, bool py, bool pz);
		void setGridDimensions(Coordinate newGridSize);
		void setNumBins(int newNumBinX, int newNumBinY, int newNumBinZ);
        void setBinSize(double desiredBinSize);
		void setPeriodic(bool newPeriodicX, bool newPeriodicY, bool newPeriodicZ);
		int putInGrid(int tag, Coordinate c);
		bool removeFromGrid(int tag);
		void clearGrid();
        std::tuple<int, int, int> binToGridPositions(int binIndex);
        int gridPositionsToBin(int cellX, int cellY, int cellZ);
	
		std::vector<int> getNeighboringTags(int particleTag);
        std::vector<int> getNeighboringTags(int particleTagI, int particleTagJ);
		int getNumBins();
        
        std::tuple<int, int, int, int> getBin(Coordinate c);
        
        std::vector<int> getBinsInRegion(Coordinate centerPosition, double radius);
        std::vector<int> getTagInBins(std::vector <int> binList);
        std::vector<int> getTagsInRegion(Coordinate centerPosition, double radius);
        
        std::vector<int> getTagsAboveZVal(double zPos);
        std::vector<int> getTagsBetweenZVals(double zPosLow, double zPosHigh);
};

#endif