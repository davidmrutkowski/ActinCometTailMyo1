/**
 * Copyright (C) 2024 Lehigh University.
 *
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see http://www.gnu.org/licenses/.
 *
 * Author: David Rutkowski (dmr518@lehigh.edu)
 */

#include "Simulation.h"

Simulation::Simulation()
{
    // default values
    temperature = 300.0;
    
    nucleatorBeadRadius = 1.5*0.5;
    
    eta = 0.301;
    numMonomersPerBead = 0;
    dt = 0.1;
    totalSimulationTime = 1.0;
    snapshotTime = 0.01;
    persistenceLength = 17.0;
    
    nucleatorFilamentInteractionDist = rRepulsiveInteraction*0.5 + nucleatorBeadRadius;
    
    numThreads = 1;
    thermalForcesOn = true;
    branchAngleK = 0.1;
    filamentBondK = 100.0;
    pullForceMag = 0.00135*0.5*3;
    // diameter of rod
    rRepulsiveInteraction = 20e-3;
    severCrosslinksWithExtension = false;
    crosslinkNeighboringUponAdding = true;
    filamentBondsInBondList = true;
    excludedVolOn = false;
    backPullingForceOn = false;
    uniformPullingForceOn = false;
    
    FABeadPos = {0.0, 0.0, -1.0-0.125-0.5};
    FABeadRadius = 0.25*0.5;
    FABeadLength = 1.25 - FABeadRadius * 2;
    
    nascentFABindProbability = 0.5;
    nascentFAUnbindProbability = 0.5;
    
    severingDistanceMultiplier = 10.0;
    
    forceMyosin = 0.0;
    
    
    branchingRatePerSegment = 3.5*0.5*1.0*0.3 *(40.0/140.0);
    
    
    currSimTime = 0.0;
    
    outputFileName = "";
    
    initDerivedValues();

    //https://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution
    seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::cout << "Seed used: " << seed << std::endl;
    dis = std::uniform_real_distribution<> (0.0, 1.0);
    gaussianDis = std::normal_distribution<> (0.0, 1.0);
    
    //https://stackoverflow.com/questions/37305104/how-to-generate-random-numbers-in-a-thread-safe-way-with-openmp
	for (int i = 0; i < omp_get_max_threads(); i++)
	{
        std::mt19937 gen(seed+i);
		generators.emplace_back(std::mt19937(gen));
	}
}

void Simulation::initDerivedValues()
{
    L_zero = 2.7e-3 * numMonomersPerBead;
    zeta = 4*pi*eta*L_zero / (0.84 + log(L_zero/2/0.0035));
    snapshotStep = (int)round(snapshotTime / dt);
    numSteps = (int)round(totalSimulationTime / dt);
    filamentK = temperature * boltzmannConstant * persistenceLength;
    thermal_force = sqrt(2*boltzmannConstant*temperature*zeta/dt);
    zetaOfNucleatorBead = 6.0 * pi * eta * nucleatorBeadRadius * 200.0;
    
    std::cout << zetaOfNucleatorBead << std::endl;
    
    thermal_force_inNascentFA = sqrt(2*boltzmannConstant*temperature*zetaOfNucleatorBead/dt);
    delXDividedByKT = delX / boltzmannConstant / temperature;
}

void Simulation::readParameterFile(std::string fileName)
{
    ifstream inputfile (fileName);
    double timeBetweenDeltaSteps = 0.0;
    
    if (inputfile.is_open())
	{
		int linecount = 0;
		
        std::string line;
              
		while(std::getline(inputfile, line))
		{	
			std::istringstream iss(line);
			
			std::string category;
            std::string value;

            iss >> category >> value;

            if(category.compare("dt") == 0)
            {
                dt = std::stod(value);
            }
            else if(category.compare("totalSimulationTime") == 0)
            {
                totalSimulationTime = std::stod(value);
            }
            else if(category.compare("snapshotTime") == 0)
            {
                snapshotTime = std::stod(value);
            }
            else if(category.compare("persistenceLength") == 0)
            {
                persistenceLength = std::stod(value);
            }
            else if(category.compare("numMonomersPerBead") == 0)
            {
                numMonomersPerBead = std::stod(value);
            }
            else if(category.compare("temperature") == 0)
            {
                temperature = std::stod(value);
            }
            else if(category.compare("eta") == 0)
            {
                eta = std::stod(value);
            }
            else if(category.compare("positions") == 0)
            {
                positionFileName = value;
            }
            else if(category.compare("bonds") == 0)
            {
                bondFileName = value;
            }
            else if(category.compare("angles") == 0)
            {
                angleFileName = value;
            }
            else if(category.compare("numThreads") == 0)
            {
                numThreads = std::stod(value);
            }
            else if(category.compare("thermalForcesOn") == 0)
            {
                istringstream(value) >> std::boolalpha >> thermalForcesOn;
            }
            else if(category.compare("filamentBondK") == 0)
            {
                filamentBondK = std::stod(value);
            }
            else if(category.compare("branchAngleK") == 0)
            {
                branchAngleK = std::stod(value);
            }
            else if(category.compare("severCrosslinksWithExtension") == 0)
            {
                istringstream(value) >> std::boolalpha >> severCrosslinksWithExtension;
            }
            else if(category.compare("outputName") == 0)
            {
                outputFileName = value;
            }
            else if(category.compare("excludedVolOn") == 0)
            {
                istringstream(value) >> std::boolalpha >> excludedVolOn;
            }
            else if(category.compare("backPullingForceOn") == 0)
            {
                istringstream(value) >> std::boolalpha >> backPullingForceOn;
            }
            else if(category.compare("uniformPullingForceOn") == 0)
            {
                istringstream(value) >> std::boolalpha >> uniformPullingForceOn;
            }
            else if(category.compare("etaInNascentAdhesion") == 0)
            {
                etaInNascentAdhesion = std::stod(value);
            }
            else if(category.compare("nascentFABindProbability") == 0)
            {
                nascentFABindProbability = std::stoi(value);
            }
            else if(category.compare("nascentFAUnbindProbability") == 0)
            {
                nascentFAUnbindProbability = std::stoi(value);
            }
            else if(category.compare("reactionRateNascentFABind") == 0)
            {
                double reactionRate = std::stoi(value);
                nascentFABindProbability = 1.0 - exp(-reactionRate * dt);
            }
            else if(category.compare("reactionRateNascentFAUnbind") == 0)
            {
                double reactionRate = std::stoi(value);
                nascentFAUnbindProbability = 1.0 - exp(-reactionRate * dt);
            }
            else if(category.compare("maxLength") == 0)
            {
                maxLength = std::stoi(value);
            }
            else if(category.compare("forceMyosin") == 0)
            {
                forceMyosin = std::stod(value);
            }
            else
            {
                throw std::runtime_error("Unknown value in input file associated with: " + category);
            }
            
            linecount++;
			
		}
        
		inputfile.close();
    }
    else
    {
        throw std::runtime_error("main could not open: " + fileName);
    }
    
    initDerivedValues();
    
    Coordinate boxSize = {20.0, 20.0, 20.0};
    bool periodicX = false;
    bool periodicY = false;
    bool periodicZ = false;
    
    simBox = SimulationBox(boxSize, periodicX, periodicY, periodicZ);
    particleGrid = Grid(boxSize, periodicX, periodicY, periodicZ);
    particleGrid.setBinSize(L_zero*0.5);
    
    if(!positionFileName.empty())
    {
        readXYZ();
    }
    
    if(bondFileName.empty() && !positionFileName.empty())
    {
        
        // assumes that bond and angle file names are the same as positionFileName just with a different extension
        bondFileName = positionFileName.substr(0, positionFileName.length() - 4) + ".bnd";
        
        cout << "Assuming bond file is named: " << bondFileName << endl;
    }
    readBondFile();
    
    cout << "read bond file" << endl;
    
    if(angleFileName.empty() && !positionFileName.empty())
    {
        angleFileName = positionFileName.substr(0, positionFileName.length() - 4) + ".ang";
        
        cout << "Assuming angle file is named: " << angleFileName << endl;
    }
    readAngleFile();
    
    double updateGridTime = 0.05*(140.0/40.0);
    updateGridStep = (int)round(updateGridTime / dt);
    
    criticalBendingForceMag = 10.0;
    criticalBondForceMag = -50.0;
    
    minCrosslinkBondDist = 0.030;
    maxCrosslinkBondDist = 0.040;
    
    // filament bonds, type 0
    bondTypes.emplace_back(SpringBondInfo(filamentBondK, L_zero));
    
    // permanent crosslink bonds, type 1
    bondTypes.emplace_back(SpringBondInfo(100.0, 0.5*(minCrosslinkBondDist + maxCrosslinkBondDist)));
    
    // temporary crosslink bonds, type 2
    bondTypes.emplace_back(SpringBondInfo(100.0, 0.5*(minCrosslinkBondDist + maxCrosslinkBondDist)));
    
    // filament bonds, for arp2/3, type 3
    bondTypes.emplace_back(SpringBondInfo(filamentBondK, rRepulsiveInteraction / sin(70.0*pi/180.0)));
    
    // growing segment lengths, additional types, types [4 to 4+numMonomersPerBead-2]
    for(int i = 1; i < numMonomersPerBead; i++)
    {
        double tmpBondK = filamentBondK;        
        bondTypes.emplace_back(SpringBondInfo(tmpBondK, 2.7E-3*i));
    }
    
    thermal_force = sqrt(2*boltzmannConstant*temperature*zeta/dt);
}

void Simulation::readXYZ()
{
    ifstream inputfile (positionFileName);
    
    int numParticles = 0;
    
    if (inputfile.is_open())
	{
		int linecount = 0;
        
        int currFilamentType = -1;
        int newFilamentTag = -1;
        
        std::string line;
        
		while(std::getline(inputfile, line))
		{
            if(line.size())
            {
                std::istringstream iss(line);
                
                if(linecount == 0)
                {
                    iss >> numParticles;
                }
                else if(linecount == 1)
                {
                    std::string tempString;
                    iss >> tempString;
                    
                    int pos = tempString.find("=");
                    
                    std::string timeString = tempString.substr(pos+1,tempString.length());
                    
                    currSimTime = std::stod(timeString);
                }
                else if(linecount > 1)
                {
                    int tag, type;
                    double x,y,z;
                    int viscosityType;
                    int cappedStatus;
                    
                    iss >> type >> tag >> x >> y >> z >> cappedStatus;
                    
                    type = type - 1;
                    
                    
                    int filamentIndex = f_TaggedVector.getIndexOfTag(type);
                    
                    
                    if(filamentIndex < 0 && type >= 0)
                    {
						if(filaments.size() > 0 && filaments[filaments.size()-1].getNumParticles() <= 2)
						{
							std::cout << "this filament has less than 3 segments " << filaments[filaments.size()-1].getNumParticles() << " " << type << std::endl;
						}
						
                        Filament newFil;
                        newFil.setIsDaughter(false);
                        filaments.push_back(newFil);
                        newFilamentTag = f_TaggedVector.add();
                        
                        if(newFilamentTag != type)
                        {
                            f_TaggedVector.setTagAtPos(filaments.size()-1, type);
                        }
                        
                        filamentIndex = filaments.size() - 1;
                    }

                    Coordinate newCoordinate = {x, y, z};
                    
                    Particle newParticle(x, y, z);
                    newParticle.setType(type);
                    
                    int numStressMeasurements = 0;
                    
                    int newTag = pinfo.addParticle(newParticle, numStressMeasurements);
                    if(newTag != tag)
                    {
                        // set tag of this particle to its correct value listed in the xyz file
                        pinfo.setTagAtPos(pinfo.getNumParticles()-1, tag);
                        newTag = tag;
                    }
                    
                    if(type >= 0)
                    {
                        if(cappedStatus == 1)
                        {
                            //need to figure out if this bead is the pointed or barbed end
                            if(filaments[filaments.size()-1].getNumParticles() <= 1)
                            {
                                // assume this is the pointed end (first bead in the filament)
                                // this shouldn't happen currently since pointed ends are not capped currently
                                filaments[filaments.size()-1].setCappedPointedEnd(true);
								std::cout << "pointed end capped on restart!" << std::endl;
                            }
                            else
                            {
                                filaments[filaments.size()-1].setCappedBarbedEnd(true);
                            }
                        }
                        
                        int priorTag = filaments[filamentIndex].addParticleBack(newTag);
                        
                        if(!filamentBondsInBondList && priorTag > -1)
                        {
                            // adds a bond between priorTag and newTag of type 0 if they are on the same filament
                            // as indicated from the type column in the xyz file
                            Bond newBond = {priorTag, newTag, 0};
                            
                            pinfo.addBond(newBond);
                        }
                    }
                }
                
                linecount++;
            }
        }
        
        if(linecount-2 != numParticles)
        {
            // warning that the number of particles in the simulation is different from numParticles
            std::cout << "Error. " << numParticles << " particles desired in simulation but " << linecount-2 << " particles in " << this->positionFileName;
        }
        
        inputfile.close();
    }
    else
    {
        cout << "Could not open " + positionFileName << endl;
    }
}

void Simulation::readBondFile()
{
    ifstream inputfile (bondFileName);
    if (inputfile.is_open())
	{
		int linecount = 0;
		
        int numBonds = 0;
        int countBonds = 0;
        
        int numParticles = pinfo.getNumParticles();
        
        std::string line;
        
        int maxCurrTag = pinfo.getCurrMaxTag();
        
		while(std::getline(inputfile, line))
		{
            if(line.size())
            {
                std::istringstream iss(line);
                
                if(linecount == 0)
                {
                    iss >> numBonds;
                }
                else if(linecount > 1)
                {
                    int particleI, particleJ, bondType;
                    double bond_creationTime, bond_destructionTime;
                    
                    iss >> bondType >> particleI >> particleJ >> bond_creationTime >> bond_destructionTime;
                    
                    if(particleI < 0 || particleI > maxCurrTag || particleJ < 0 || particleJ > maxCurrTag)
                    {
                        // warning about particleI and/or particleJ being out of bounds
                        throw std::runtime_error(std::to_string(particleI) + " or " + std::to_string(particleJ) + " tags are larger than current max tag " + std::to_string(maxCurrTag));
                    }
                    
                    //add bond tag to both particleI and particleJ
                    pinfo.addBond(Bond {particleI, particleJ, bondType, bond_creationTime, bond_destructionTime});
                    countBonds++;
                }
                
                linecount++;
            }
        }

        if(countBonds != numBonds)
        {
            // warning that the number of particles in the simulation is different from numBonds
            std::cout << "Error. " << numBonds << " bonds desired in simulation but " << countBonds << " bonds in " << this->bondFileName;
        }
        
        inputfile.close();
    }
}

void Simulation::readAngleFile()
{
    ifstream inputfile (angleFileName);
    if (inputfile.is_open())
	{
		int linecount = 0;
		
        int numAngles = 0;
        int countAngles = 0;
        
        std::string line;
        
		while(std::getline(inputfile, line))
		{
            if(line.size())
            {
                std::istringstream iss(line);
                
                if(linecount == 0)
                {
                    iss >> numAngles;
                }
                else if(linecount > 1)
                {
                    int particleI, particleJ, particleK, angleType;
                    iss >> angleType >> particleI >> particleJ >> particleK;
                    
                    pinfo.addAngle(Angle {particleI, particleJ, particleK, angleType});
                    countAngles++;
                }
                
                linecount++;
            }
        }
        
        if(countAngles != numAngles)
        {
            // warning that the number of angles in the simulation is different from numAngles
            throw std::runtime_error("Number of angles desired in simulation is " + std::to_string(numAngles) + " but " + std::to_string(countAngles) + " angles in " + angleFileName);
        }
        
        inputfile.close();
    }
}

void Simulation::moveParticles()
{
    int numParticles = pinfo.getNumParticles();
    
    #pragma omp parallel for
    for(int i = 0; i < numParticles; i++)
    {
        Coordinate currPosition = pinfo.getPosByIndex(i);
        Coordinate currForce = pinfo.getForceByIndex(i);
        
        double vx, vy, vz;
        
        if(i == 0)
        {
            // is nucleator bead
            vx = currForce.x / zetaOfNucleatorBead;
            vy = currForce.y / zetaOfNucleatorBead;
            vz = currForce.z / zetaOfNucleatorBead;
        }
        else
        {
            // not nucleator bead
            vx = currForce.x / zeta;
            vy = currForce.y / zeta;
            vz = currForce.z / zeta;
        }
    
                
        Coordinate nextPosition;
        nextPosition.x = currPosition.x + vx * dt;
        nextPosition.y = currPosition.y + vy * dt;
        nextPosition.z = currPosition.z + vz * dt;
        
        if(std::isnan(nextPosition.x))
            throw std::runtime_error("Particle " + std::to_string(i) + " has a NaN position in the x direction");
        
        nextPosition = simBox.periodicWrap(nextPosition);
        
        pinfo.setPosByIndex(i, nextPosition);
	}
}

void Simulation::putAllParticlesInGrid()
{
    particleGrid.clearGrid();
    
    int numParticles = pinfo.getNumParticles();
    for(int i = numParticles-1; i > -1; i--)
    {
        int tempTag = pinfo.getTagAtIndex(i);
        
        // if this tag exists
        if(tempTag >= 0)
        {
            int newBin = particleGrid.putInGrid(tempTag, pinfo.getPosByIndex(i));
        }
    }
    // set new particle neighbors
    #pragma omp parallel for
    for(int p = 0; p < pinfo.getNumParticles(); p++)
    {
        pinfo.setNeighborsByIndex(p, particleGrid.getNeighboringTags(pinfo.getTagAtIndex(p)));
    }
    
    std::vector <int> possibleBondNeighborTypes = {0,3};
    // now set bond neighbors but only for bonds with type 0 (filament bonds are only included for excluded volume)
    if(excludedVolOn)
    {
        #pragma omp parallel for
        for(int b = 0; b < pinfo.getNumBonds(); b++)
        {
            pinfo.resetBondNeighborTagsByIndex(b, possibleBondNeighborTypes);
            
            // calc bond distance here between b and neighbors and if they are greater than some value then remove them
            Bond bondB = pinfo.getBondByIndex(b);
            Coordinate b_i = pinfo.getPos(bondB.i);
            Coordinate b_j = pinfo.getPos(bondB.j);
                
            std::vector <int> newNeighboringTags = pinfo.getBondNeighborTagsByIndex(b);
            
            for(int neigh = newNeighboringTags.size()-1; neigh > -1; neigh--)
            {
                Bond neighborBond = pinfo.getBond(newNeighboringTags[neigh]);
                Coordinate neigh_i = pinfo.getPos(neighborBond.i);
                Coordinate neigh_j = pinfo.getPos(neighborBond.j);
                
                Coordinate dummy_k, dummy_l;

                Coordinate d_kl = closestDistanceBetweenLines(b_i, b_j, neigh_i, neigh_j, dummy_k, dummy_l);
                double currMinDist = d_kl.getMagnitude();
                if(currMinDist > rRepulsiveInteraction*3.0)
                {
                    // remove this bond from newNeighboringTags
                    newNeighboringTags[neigh] = newNeighboringTags.back();
                    newNeighboringTags.pop_back();
                }
            }

            
            pinfo.setBondNeighborTagsByIndex(b, newNeighboringTags);
        }
    }
}

void Simulation::calcFilamentForces(int f)
{
    std::deque <int> orgTagList = filaments[f].getAllTags();
    
    int numTags = orgTagList.size();
    int removedTags = 0;
    
    Coordinate ab, bc;
    Coordinate abUnit, bcUnit;
    int aTag, bTag, cTag;
    Coordinate posBTag, posCTag;
    double abMag, bcMag;
    
    std::mt19937& threadEngine = generators[omp_get_thread_num()]; 
        
        
    
    for(int p = 0; p < numTags-1; p++)
    {
        // shift by one bead to the left since last loop calculated these parameters
        ab = bc;
        aTag = bTag;
        abUnit = bcUnit;
        abMag = bcMag;
        
        if(p == 0)
        {
            bTag = orgTagList[p];
            posBTag = pinfo.getPos(bTag);
        }
        else
        {
            bTag = cTag;
            posBTag = posCTag;
        }

        cTag = orgTagList[p+1];
        
        posCTag = pinfo.getPos(cTag);
        
        bc = simBox.calcDisplacement(posBTag, posCTag);
        bcMag = bc.getMagnitude();
        bcUnit = bc / bcMag;
        
        // assumes that all bonds in a filament are the same and 
        // the parameters for this bond are located in the first position of vector bondTypes
        
        // search for bond tag connecting b & c
        
        int tempBondTag = pinfo.findBondTag(bTag, cTag);
        
        Bond tempBond = pinfo.getBond(tempBondTag);
        
        
        double bondForce = bondTypes[tempBond.bondTypeIndex].calcForce(bcMag);
        
        bool brokeDueToBending = false;
        if(p > 0)
        {   
            // calc bending forces
            // dont want to do this if just deleted the previous bond
            Coordinate cb = -bc;
            
            double averageMag = (abMag + bcMag) * 0.5;
            double invAverageMag = 1.0 / averageMag;
            
            double dotProd = ab.x*cb.x + ab.y*cb.y + ab.z*cb.z;

            Coordinate forceA = -filamentK * invAverageMag * (cb / bcMag / abMag + dotProd / bcMag * (-ab / (abMag*abMag*abMag)));
            Coordinate forceC = -filamentK * invAverageMag * (ab / abMag / bcMag + dotProd / abMag * (-cb / (bcMag*bcMag*bcMag)));
            
            
            Coordinate forceB = forceA + forceC;
            
            double tempDot = -abUnit.x*bcUnit.x - abUnit.y*bcUnit.y - abUnit.z*bcUnit.z;
                    
            if(tempDot < -1.0)
            {
                tempDot = -1.0;
            }
            else if(tempDot > 1.0)
            {
                tempDot = 1.0;
            }
            
            // if forceB mag is greater than some value break filament
            double theta = acos(tempDot);
            
            {
                pinfo.addForce(aTag, forceA);
                pinfo.addForce(bTag, -forceB);
                pinfo.addForce(cTag, forceC);
            }
        }
        
        Coordinate tempForce = bondForce * bcUnit;
        
        double currBondDestructionTime = 0.0;
        if(tempBondTag >= 0)
        {
            currBondDestructionTime = pinfo.getBond(tempBondTag).destructionTime;
        }
        else
        {            
            throw std::runtime_error("Could not find bondTag between particles with tags: " + std::to_string(bTag) + " " + std::to_string(cTag));
        }

        
        if(tempBond.bondTypeIndex == 0 && (currSimTime >= currBondDestructionTime || bondForce < criticalBondForceMag && brokeDueToBending == false))
        {
            //instead add this tag to the flaggedBondList
            
            #pragma omp critical
            {
				flaggedBondList.push_back(tempBondTag);
            }
        }
        else if(brokeDueToBending == false)
        {            
            pinfo.addForce(bTag, tempForce);
            pinfo.addForce(cTag, -tempForce);
        }
    }    
}

double Simulation::calcMinDistanceCylinders(Coordinate ri, Coordinate rj, Coordinate ri2, Coordinate rj2)
{
    Coordinate R_k = simBox.calcDisplacement(ri, rj);
    
    Coordinate P_k = 0.5 * R_k + rj;
    P_k = simBox.periodicWrap(P_k);
    
    double R_kSquared = R_k*R_k;

    Coordinate R_l = simBox.calcDisplacement(ri2, rj2);

    Coordinate P_l = simBox.periodicWrap(0.5 * R_l + rj2);
    
    double R_lSquared = R_l*R_l;
    
    double R_klSquared = R_k*R_l;
    
    Coordinate firstTerm = simBox.calcDisplacement(P_k, P_l);
    
    Coordinate secondTerm = simBox.periodicWrap(R_lSquared * R_k - R_klSquared*R_l);
    
    double tk = (firstTerm) * (secondTerm);
    tk = tk / (R_klSquared*R_klSquared - R_kSquared * R_lSquared);
    
    if(tk > 0.5)
        tk = 0.5;
    else if(tk < -0.5)
        tk = -0.5;
    
    firstTerm = -firstTerm;
    
    secondTerm = simBox.periodicWrap(R_kSquared * R_l - R_klSquared*R_k);
    
    double tl = (firstTerm) * (secondTerm);
    tl = tl / (R_klSquared*R_klSquared - R_kSquared * R_lSquared);
    
    if(tl > 0.5)
        tl = 0.5;
    else if(tl < -0.5)
        tl = -0.5;
    
    Coordinate closestPointOnK = simBox.periodicWrap(P_k + tk * R_k);
    
    Coordinate closestPointOnL = simBox.periodicWrap(P_l + tl * R_l);

    // points towards k from l rather than the reverse
    Coordinate d_kl = simBox.calcDisplacement(closestPointOnK, closestPointOnL);
    
    double minDist = d_kl.getMagnitude();					
    
    return minDist;
}

double Simulation::calcForces()
{
    int numParticles = pinfo.getNumParticles();

    // make tempForces vectors larger if needed and set new elements to zero
    // tempStress currently only holds one type of stress (bond)!
    if(tempForces[0].size() != numParticles)
    {
        //#pragma omp parallel for
        for(int thread = 0; thread < omp_get_max_threads(); thread++)
        {
            int oldSize = tempForces[thread].size();
            
            tempForces[thread].resize(pinfo.getNumParticles());
            for(int i = oldSize; i < pinfo.getNumParticles(); i++)
            {
                tempForces[thread][i] = Coordinate {0.0, 0.0, 0.0};
            }
        }
    }
    
    #pragma omp parallel for
    for(int p = 0; p < numParticles; p++)
    {
        Coordinate currPosition = pinfo.getPosByIndex(p);

        std::mt19937& threadEngine = generators[omp_get_thread_num()];
        
        int currTag = pinfo.getTagAtIndex(p);
        std::vector <int> iNeighbors = pinfo.getNeighbors(currTag);
        
        // loop over extra bonded particles
        std::vector <struct HalfBond> tempBondList = pinfo.getBondTagsByIndex(p);
        for(int b = 0; b < tempBondList.size(); b++)
        {
            Bond firstBond = pinfo.getBond(tempBondList[b].bondTag);
            
            // this was wrong in ASCB version
            if(currTag < tempBondList[b].j)
            {
                // bondTag is  being used here as bondType and bond list is not used at all in pinfo
                
                Coordinate ri = pinfo.getPos(currTag);
                Coordinate rj = pinfo.getPos(firstBond.j);
                Coordinate rij = simBox.calcDisplacement(ri, rj);
                
                
                // excluded volume with nucleator bead
                Coordinate closestK;
                    
                Coordinate d_kl = closestDistanceBetweenLineAndPoint(ri, rj, pinfo.getPos(0), closestK);
                
                double minDist = d_kl.getMagnitude();					

                if(minDist < nucleatorFilamentInteractionDist)
                {
                    double forceMag = 0.0;
                    
                    double kr = -1000.0;

                    forceMag = kr * (minDist - nucleatorFilamentInteractionDist);

                    Coordinate Fr = forceMag * d_kl.getUnitCoord();
                    
                    double length_k = rij.getMagnitude();
                    
                    double firstSegmentLength = simBox.calcDistance(closestK, rj);
                    double secondSegmentLength = length_k - firstSegmentLength;
                    
                    Coordinate Falpha = firstSegmentLength / length_k * (Fr);
                    Coordinate Fbeta = secondSegmentLength / length_k * (Fr);
                    
                    tempForces[omp_get_thread_num()][pinfo.getIndexOfTag(currTag)] = tempForces[omp_get_thread_num()][pinfo.getIndexOfTag(currTag)] + Falpha;
                    tempForces[omp_get_thread_num()][pinfo.getIndexOfTag(tempBondList[b].j)] = tempForces[omp_get_thread_num()][pinfo.getIndexOfTag(tempBondList[b].j)] + Fbeta;
                    
                    tempForces[omp_get_thread_num()][0] = tempForces[omp_get_thread_num()][0] - Fr;
                }
                // end of excluded volume
                
                
                // if this bond is a crosslink bond then calculate bond forces on it (but no excluded forces)
                if(firstBond.bondTypeIndex > 0)
                {
                    double rijMag = rij.getMagnitude();
                    
                    double bondDestructionTime = firstBond.destructionTime;
                    
                    Coordinate rijUnit = rij / rijMag;
                    
                    double bondForce = bondTypes[firstBond.bondTypeIndex].calcForce(rijMag);
                    Coordinate bondForceVector = bondForce * rijUnit;
                
                    int currPosJ = pinfo.getIndexOfTag(tempBondList[b].j);
                
                
                    if(firstBond.bondTypeIndex == 3 && (currSimTime >= firstBond.destructionTime || bondForce < -30.0))
                    {
                        //instead add this tag to the flaggedBondList
                        
                        #pragma omp critical
                        {
							flaggedBondList.push_back(tempBondList[b].bondTag);
                        }
                    }
        
                    tempForces[omp_get_thread_num()][p] = tempForces[omp_get_thread_num()][p] + bondForceVector;
                    tempForces[omp_get_thread_num()][currPosJ] = tempForces[omp_get_thread_num()][currPosJ] - bondForceVector;
                }
                
            }
        }
        
        // loop over extra angle particles
        std::vector <struct Angle> tempAngleList = pinfo.getAngleTagsByIndex(p);
        for(int la = tempAngleList.size()-1; la > -1; la--)
        {
            if(currTag == tempAngleList[la].i)
            {
                Angle ang = pinfo.getAngle(tempAngleList[la].angleTypeIndex);
                
                if(ang.angleTypeIndex == 0)
                {
                    // 70 degree angle
                    
                    double theta0 = cos(1.22173);
                    double angleK = branchAngleK;
                    
                    Coordinate a = pinfo.getPos(ang.i);
                    Coordinate b = pinfo.getPos(ang.j);
                    Coordinate c = pinfo.getPos(ang.k);
                    
                    Coordinate ab = simBox.calcDisplacement(a, b);
                    double abMag = ab.getMagnitude();
                    Coordinate abUnit = ab / abMag;
                    
                    Coordinate cb = simBox.calcDisplacement(c, b);
                    double cbMag = cb.getMagnitude();
                    Coordinate cbUnit = cb / cbMag;
                    
                    double tempDot = abUnit.x*cbUnit.x + abUnit.y*cbUnit.y + abUnit.z*cbUnit.z;
                    
                    if(tempDot < -1.0)
                    {
                        tempDot = -1.0;
                    }
                    else if(tempDot > 1.0)
                    {
                        tempDot = 1.0;
                    }
                    
                    double theta = acos(tempDot);
                    
                    double dotProd = ab.x*cb.x + ab.y*cb.y + ab.z*cb.z;
                    
                    Coordinate forceA = -2.0*angleK * (tempDot - theta0) * (cb / cbMag / abMag + dotProd / cbMag * (-ab / (abMag*abMag*abMag)));
                    Coordinate forceC = -2.0*angleK * (tempDot - theta0) * (ab / abMag / cbMag + dotProd / abMag * (-cb / (cbMag*cbMag*cbMag)));
                    
                    Coordinate forceB = forceA + forceC;
                    
                    int aIndex = pinfo.getIndexOfTag(ang.i);
                    int bIndex = pinfo.getIndexOfTag(ang.j);
                    int cIndex = pinfo.getIndexOfTag(ang.k);
                    
                    tempForces[omp_get_thread_num()][aIndex] = tempForces[omp_get_thread_num()][aIndex] + forceA;
                    tempForces[omp_get_thread_num()][bIndex] = tempForces[omp_get_thread_num()][bIndex] - forceB;
                    tempForces[omp_get_thread_num()][cIndex] = tempForces[omp_get_thread_num()][cIndex] + forceC;
                    
                    // angle criteria +- 25 degrees  
                    if(theta <= 0.785398 || theta >= 1.65806)
                    {
                        // then remove this branch (angles and bond)
                        
                        // both angles will be removed when flagged bond is removed
                            
                        // and remove bond between either a-b or b-c, two different filaments
                        int aFilamentType = pinfo.getBeadToFilamentByIndex(aIndex);
                        int bFilamentType = pinfo.getBeadToFilamentByIndex(bIndex);
                        int cFilamentType = pinfo.getBeadToFilamentByIndex(cIndex);
                        
                        // a and b are of the same filament, this means c is on another filament (daughter)
                        if(aFilamentType == bFilamentType)
                        {
                            //then remove b-c bond
                            std::vector <struct HalfBond> tmpBranchBondList = pinfo.getBondTagsByIndex(bIndex);
                            for(int b = 0; b < tmpBranchBondList.size(); b++)
                            {
                                Bond tmpBond = pinfo.getBond(tmpBranchBondList[b].bondTag);
                                
                                if(tmpBranchBondList[b].j == ang.k)
                                {
                                    // remove this bond
                                    #pragma omp critical
                                    {
                                        cout << "remove branch bond: " << ang.j << " " << tmpBranchBondList[b].j << " " << theta << " " << dotProd << endl;
                                        flaggedBondList.push_back(tmpBranchBondList[b].bondTag);
                                    }
                                    break;
                                }
                            }
                        }
                        else
                        {
                            // then remove a-c bond
                            std::vector <struct HalfBond> tmpBranchBondList = pinfo.getBondTagsByIndex(bIndex);
                            for(int b = 0; b < tmpBranchBondList.size(); b++)
                            {
                                Bond tmpBond = pinfo.getBond(tmpBranchBondList[b].bondTag);
                                
                                if(tmpBranchBondList[b].j == ang.i)
                                {
                                    // remove this bond
                                    #pragma omp critical
                                    {
                                        cout << "remove (else) bond: " << tmpBranchBondList[b].bondTag << endl;
                                        flaggedBondList.push_back(tmpBranchBondList[b].bondTag);
                                    }
                                    break;
                                }
                            }
                        }                        
                    }
                        
                }
                else if(ang.angleTypeIndex == 1)
                {
                    // straight angle
                    
                    Coordinate a = pinfo.getPos(ang.i);
                    Coordinate b = pinfo.getPos(ang.j);
                    Coordinate c = pinfo.getPos(ang.k);
                    
                    Coordinate ab = simBox.calcDisplacement(a, b);
                    double abMag = ab.getMagnitude();
                    Coordinate abUnit = ab / abMag;
                    
                    Coordinate cb = simBox.calcDisplacement(c, b);
                    double cbMag = cb.getMagnitude();
                    Coordinate cbUnit = cb / cbMag;
                    
                    double invAverageMag = 1.0 / 0.1;
                    
                    double dotProd = ab.x*cb.x + ab.y*cb.y + ab.z*cb.z;
                    
                    Coordinate forceA = -filamentK * invAverageMag * (cb / cbMag / abMag + dotProd / cbMag * (-ab / (abMag*abMag*abMag)));
                    Coordinate forceC = -filamentK * invAverageMag * (ab / abMag / cbMag + dotProd / abMag * (-cb / (cbMag*cbMag*cbMag)));
                    
                    Coordinate forceB = forceA + forceC;
                    
                    int aIndex = pinfo.getIndexOfTag(ang.i);
                    int bIndex = pinfo.getIndexOfTag(ang.j);
                    int cIndex = pinfo.getIndexOfTag(ang.k);
                    
                    tempForces[omp_get_thread_num()][aIndex] = tempForces[omp_get_thread_num()][aIndex] + forceA;
                    tempForces[omp_get_thread_num()][bIndex] = tempForces[omp_get_thread_num()][bIndex] - forceB;
                    tempForces[omp_get_thread_num()][cIndex] = tempForces[omp_get_thread_num()][cIndex] + forceC;
                }
            }
        }        
	}
    
    if(excludedVolOn)
    {
        int currNumBonds = pinfo.getNumBonds();
        
        std::vector <int> countPerThread (omp_get_max_threads());
        
        
        // bonds of lower indexes have larger numbers of neighbors on average so the work per thread will not be balanced
        #pragma omp parallel for schedule(static, 1)
        for(int b = 0; b < currNumBonds; b++)
        {            
            std::vector <int> currNeighborBondTags = pinfo.getBondNeighborTagsByIndex(b);
            Bond firstBond = pinfo.getBondByIndex(b);
            
            int currPosI, currPosJ;
            Coordinate ri, rj;
            
            try
            {
                currPosI = pinfo.getIndexOfTag(firstBond.i);
                currPosJ = pinfo.getIndexOfTag(firstBond.j);
                
                ri = pinfo.getPosByIndex(currPosI);
                rj = pinfo.getPosByIndex(currPosJ);
            }
            catch (const std::runtime_error& error)
            {
                #pragma omp critical
                {
                    cout << b << " " << firstBond.i << " " << firstBond.j << " " << currNumBonds << " " << currPosI << " " << currPosJ << endl;
                    exit(0);
                }
            }
                
            Coordinate R_k = simBox.calcDisplacement(ri, rj);
            double R_kSquared = R_k*R_k;
            
            Coordinate P_k = 0.5 * R_k + rj;
            P_k = simBox.periodicWrap(P_k);
            
            for(int neigh = 0; neigh < currNeighborBondTags.size(); neigh++)
            {
                // if tag of neighboring bond is less than tag at pos b then calculate force
                if(currNeighborBondTags[neigh] < pinfo.getBondTagAtPos(b))
                {
                    Bond secondBond = pinfo.getBond(currNeighborBondTags[neigh]);

                    
                    int currPosI2, currPosJ2;
                    Coordinate ri2, rj2;
                    
                    {
                    
                        try
                        {
                            currPosI2 = pinfo.getIndexOfTag(secondBond.i);
                            currPosJ2 = pinfo.getIndexOfTag(secondBond.j);
                            
                            ri2 = pinfo.getPosByIndex(currPosI2);
                            rj2 = pinfo.getPosByIndex(currPosJ2);
                        }
                        catch(const runtime_error& error)
                        {
                            #pragma omp critical
                            {
                                cout << "second" << " " << secondBond.i << " " << secondBond.j << " " << currPosI2 << " " << currPosJ2 << endl;
                                exit(0);
                            }
                        }

                        Coordinate R_l = simBox.calcDisplacement(ri2, rj2);

                        Coordinate P_l = simBox.periodicWrap(0.5 * R_l + rj2);
                        
                        double R_lSquared = R_l*R_l;
                        
                        double R_klSquared = R_k*R_l;

                        Coordinate closestK, closestL;
                        
                        Coordinate d_kl = closestDistanceBetweenLines(ri, rj, ri2, rj2, closestK, closestL);
                        double minDist = d_kl.getMagnitude();					
                        
                        if(minDist < rRepulsiveInteraction)
                        {
                            double forceMag = 0.0;
                            
                            double kr = -1690.0;
                            // excluded vol
                            forceMag = kr * (minDist - rRepulsiveInteraction);

                            Coordinate Fr = forceMag * d_kl.getUnitCoord();
                            
                            double length_k = R_k.getMagnitude();
                            
                            double firstSegmentLength = simBox.calcDistance(closestK, rj);
                            
                            double secondSegmentLength = length_k - firstSegmentLength;
                            
                            Coordinate Falpha = firstSegmentLength / length_k * (Fr);
                            Coordinate Fbeta = secondSegmentLength / length_k * (Fr);
                            
                            double length_l = R_l.getMagnitude();
                            
                            firstSegmentLength = simBox.calcDistance(closestL, rj2);
                            secondSegmentLength = length_l - firstSegmentLength;
                            
                            Coordinate Falpha2 = firstSegmentLength / length_l * (-Fr);
                            Coordinate Fbeta2 = secondSegmentLength / length_l * (-Fr);
                            
                            if(currPosI >= tempForces[omp_get_thread_num()].size() || currPosJ >= tempForces[omp_get_thread_num()].size() || currPosI2 >= tempForces[omp_get_thread_num()].size() || currPosJ2 >= tempForces[omp_get_thread_num()].size())
                            {
                                cout << "too big for tempForces: " << currPosI << " " << currPosJ << " " << currPosI2 << " " << currPosJ2 << " " << pinfo.getNumParticles() << " " << tempForces[omp_get_thread_num()].size() << endl;
                                exit(0);
                            }
                            
                            if(currPosI < 0 || currPosJ < 0 || currPosI2 < 0 || currPosJ2 < 0)
                            {
                                cout << "out of bounds of tempForces (-): " << currPosI << " " << currPosJ << " " << currPosI2 << " " << currPosJ2 << " " << pinfo.getNumParticles() << " " << tempForces[omp_get_thread_num()].size() << endl;
                                exit(0);
                            }
                            
                            tempForces[omp_get_thread_num()][currPosI] = tempForces[omp_get_thread_num()][currPosI] + Falpha;
                            tempForces[omp_get_thread_num()][currPosJ] = tempForces[omp_get_thread_num()][currPosJ] + Fbeta;
                            tempForces[omp_get_thread_num()][currPosI2] = tempForces[omp_get_thread_num()][currPosI2] + Falpha2;
                            tempForces[omp_get_thread_num()][currPosJ2] = tempForces[omp_get_thread_num()][currPosJ2] + Fbeta2;
                        }
                    }
                }
            }
        }
    }
    
    // myosin force
	std::vector<int> tagsPotentialCrosslink = particleGrid.getTagsInRegion(pinfo.getPosByIndex(0), nucleatorBeadRadius*1.25);
	
	for(int i = tagsPotentialCrosslink.size()-1; i > -1 ; i--)
	{
		double dist = simBox.calcDistance(pinfo.getPos(tagsPotentialCrosslink[i]), pinfo.getPosByIndex(0));
		
		if(((dist - nucleatorBeadRadius) > 0.025) || tagsPotentialCrosslink[i] == 0)
		{
			// remove this tag from tagsPotentialCrosslink
			tagsPotentialCrosslink[i] = tagsPotentialCrosslink.back();
			tagsPotentialCrosslink.pop_back();
		}
	}
	
	#pragma omp parallel for
	for(int p = 0; p < tagsPotentialCrosslink.size(); p++)
	{
		int currTag = tagsPotentialCrosslink[p];
		
		std::vector <struct HalfBond> tempBondList = pinfo.getBondTags(currTag);
		for(int b = 0; b < tempBondList.size(); b++)
		{
			int bond_tag_j = tempBondList[b].j;
			
			if(bond_tag_j < currTag)
			{
				Coordinate unit_vec = (simBox.calcDisplacement(pinfo.getPos(bond_tag_j), pinfo.getPos(currTag))).getUnitCoord();
				
				int currIndex = pinfo.getIndexOfTag(currTag);
				
				tempForces[omp_get_thread_num()][currIndex] = tempForces[omp_get_thread_num()][currIndex] + forceMyosin*unit_vec;
				tempForces[omp_get_thread_num()][0] = tempForces[omp_get_thread_num()][0] - forceMyosin*unit_vec;
				break;
			}
		}
	}
    //end of myosin force

    
    #pragma omp parallel for
    for(int p = 0; p < pinfo.getNumParticles(); p++)
    {
        // calculate thermal forces
        std::mt19937& threadEngine = generators[omp_get_thread_num()];
        
        Coordinate currPos = pinfo.getPosByIndex(p);
            
        if(temperature > 0.0 && thermalForcesOn == true)
        {
            // set thermal force
            if(pinfo.getViscosityTypeByIndex(p) == 0)
            {
                Coordinate c = {thermal_force*gaussianDis(threadEngine), thermal_force*gaussianDis(threadEngine), thermal_force*gaussianDis(threadEngine)};
                pinfo.setForceByIndex(p, c);
            }
            else
            {
                Coordinate c = {thermal_force_inNascentFA*gaussianDis(threadEngine), thermal_force_inNascentFA*gaussianDis(threadEngine), thermal_force_inNascentFA*gaussianDis(threadEngine)};
                pinfo.setForceByIndex(p, c);
            }
            
            
        }
        else
        {
            pinfo.setForceByIndex(p, Coordinate {0.0, 0.0, 0.0});
        }
        
        // sum up all entries in tempForces which were split across the threads
        Coordinate totalForce = {0.0, 0.0, 0.0};
        for(int thread = 0; thread < omp_get_max_threads(); thread++)
        {
            Coordinate tempForce = tempForces[thread][p];
            
            if(tempForce.x != 0.0 || tempForce.y != 0.0 || tempForce.z != 0.0)
            {
                totalForce = totalForce + tempForce;
                
                tempForces[thread][p] = Coordinate {0.0, 0.0, 0.0};
            }
        }
        pinfo.addForceByIndex(p, totalForce);
    }
    
    double totalLeadingEdgeForce = 0.0;
    // calculate filament forces (intra bond and angle forces)
    #pragma omp parallel for reduction (+:totalLeadingEdgeForce)
    for(int f = 0; f < filaments.size(); f++)
    {
        calcFilamentForces(f);
    }
    
    return 0.0;
}

//removes bonds marked for removal, returns the number of crosslinks removed
int Simulation::removeFlaggedBonds()
{  
    int numBondsRemoved = 0;
    
    while(flaggedBondList.size() > 0)
    {
        int bTag = flaggedBondList.back();
        flaggedBondList.pop_back();
        int currIndex = pinfo.getBondPosOfTag(bTag);
        
        // if this bond exists still then remove it
        if(currIndex >= 0)
        {
            Bond currBond = pinfo.getBondByIndex(currIndex);
            
            std::cout << "Removing bond with tag: " << bTag << " " << "destruction time: " << currBond.destructionTime << std::endl;
            // this bond exists
            //if(currType < 0)
            {
                //if(currType == -1)
                numBondsRemoved++;
                
                removeBondByIndex(currIndex);
            }
        }
    }
    
    return numBondsRemoved;
}

void Simulation::writeXYZFile(ofstream& outputXYZ)
{
    outputXYZ << pinfo.getNumParticles() << endl;
    outputXYZ << "t=" << currSimTime << endl;
    
    //nucleator bead
    for(int p = 0; p < 1; p++)
    {
        Coordinate tempCoordinate = pinfo.getPosByIndex(p);
            
        Coordinate tempForce = pinfo.getForceByIndex(p); 
    
        outputXYZ << 0 << " " << pinfo.getTagAtIndex(p) << " " << tempCoordinate.x << " " << tempCoordinate.y << " " << tempCoordinate.z;
        outputXYZ << " " << 0;
        

        outputXYZ << endl;
    }
    
    for(int f = 0; f < filaments.size(); f++)
    {
        int numParticlesInFilament = filaments[f].getNumParticles();
        bool isBarbedEndCapped = filaments[f].getCappedBarbedEnd();
        bool isPointedEndCapped = filaments[f].getCappedPointedEnd();
        
        for(int f2 = 0; f2 < numParticlesInFilament; f2++)
        {
            int p = pinfo.getIndexOfTag(filaments[f].getTagAtIndex(f2));
            Coordinate tempCoordinate = pinfo.getPosByIndex(p);
            
            Coordinate tempForce = pinfo.getForceByIndex(p); 
        
            outputXYZ << pinfo.getBeadToFilamentByIndex(p)+1 << " " << pinfo.getTagAtIndex(p) << " " << tempCoordinate.x << " " << tempCoordinate.y << " " << tempCoordinate.z;
            
            
            // also write out capped state of this bead
            
            if((f2 == 0 && isPointedEndCapped) || (f2 == numParticlesInFilament-1 && isBarbedEndCapped))
            {
                outputXYZ << " " << 1;
            }
            else
            {
                outputXYZ << " " << 0;
            }
            

            outputXYZ << endl;
        }
    }
}

void Simulation::writeBondFile(ofstream& outputBND)
{
    int numBonds = pinfo.getNumBonds();
    
    //if(numBonds > 0)
    {
        outputBND << numBonds << endl;
        outputBND << "t=" << currSimTime << endl;
        for(int i = 0; i < numBonds; i++)
        {
            Bond tempBond = pinfo.getBondByIndex(i);
            outputBND << tempBond.bondTypeIndex << " " << tempBond.i << " " << tempBond.j << " " << tempBond.creationTime << " " << tempBond.destructionTime << endl;
        }
        
        outputBND << endl;
    }
}

void Simulation::writeAngleFile(ofstream& outputANG)
{
    int numAngles = pinfo.getNumAngles();
    
    outputANG << numAngles << endl;
    outputANG << "t=" << currSimTime << endl;
    for(int i = 0; i < numAngles; i++)
    {
        Angle tempAngle = pinfo.getAngleByIndex(i);
        outputANG << tempAngle.angleTypeIndex << " " << tempAngle.i << " " << tempAngle.j << " " << tempAngle.k << endl;
    }
    
    outputANG << endl;
}

int Simulation::addParticlePutInGrid(Particle p, int numStressMeasurements)
{
    int newTag = pinfo.addParticle(p, numStressMeasurements);
    
    int newGrid = particleGrid.putInGrid(newTag, pinfo.getPos(newTag));
    
    std::vector <int> tempNeighbors = particleGrid.getNeighboringTags(newTag);
    for(int i = 0; i < tempNeighbors.size(); i++)
    {
        pinfo.addNeighbor(tempNeighbors[i], newTag);
    }
    
    // also need to set neighbors of newTag as well!
    pinfo.setNeighbors(newTag, tempNeighbors);
    
    return newTag;
}

int Simulation::addBondFindNeighbors(Bond newBond)
{
    int newBondTag = pinfo.addBond(newBond);
    
    Coordinate b_i = pinfo.getPos(newBond.i);
    Coordinate b_j = pinfo.getPos(newBond.j);
    
    std::vector <int> possibleBondNeighborTypes = {0,3};
    
    //sets initial bond neighbors of newly added bond
    pinfo.resetBondNeighborTagsByIndex(pinfo.getNumBonds()-1, possibleBondNeighborTypes);
    
    std::vector <int> newNeighboringTags = pinfo.getBondNeighborTagsByIndex(pinfo.getNumBonds()-1);
    
    // should measure distance between this new bond and the nucleator bead and if it is too close then remove it
    Coordinate nucleator_pos = pinfo.getPos(0);
    
    Coordinate dummy_k;

        
    Coordinate d_kl_nucleator = closestDistanceBetweenLineAndPoint(b_i, b_j, nucleator_pos, dummy_k);
                    
    if(d_kl_nucleator.getMagnitude() < nucleatorFilamentInteractionDist*0.25)
    {
        // newBond is too close to nucleator, remove this bond
        removeBond(newBondTag);
        return 1;
    }
    
    for(int neigh = newNeighboringTags.size()-1; neigh > -1; neigh--)
    {
        Bond neighborBond = pinfo.getBond(newNeighboringTags[neigh]);
		
        Coordinate neigh_i = pinfo.getPos(neighborBond.i);
        Coordinate neigh_j = pinfo.getPos(neighborBond.j);
        
        Coordinate dummy_l;

        Coordinate d_kl = closestDistanceBetweenLines(b_i, b_j, neigh_i, neigh_j, dummy_k, dummy_l);
        double currMinDist = d_kl.getMagnitude();
        if(currMinDist > rRepulsiveInteraction*5.0)
        {
            // remove this bond from newNeighboringTags
            newNeighboringTags[neigh] = newNeighboringTags.back();
            newNeighboringTags.pop_back();
        }
        else if(currMinDist < rRepulsiveInteraction*0.25)
        {
            // newBond is too close to an existing bond, remove it
            removeBond(newBondTag);
            return 2;
        }
    }
    

    pinfo.setBondNeighborTagsByIndex(pinfo.getNumBonds()-1, newNeighboringTags);
    
    // append this bond to all of newNeighboringTags
    for(int i = 0; i < newNeighboringTags.size(); i++)
    {
        int neighboringBondTag = newNeighboringTags[i];
        
        pinfo.appendBondNeighborTag(neighboringBondTag, newBondTag);
    }
    
    return 0;
}

int Simulation::addBondFindNeighborsNoRemoval(Bond newBond)
{
    int newBondTag = pinfo.addBond(newBond);
    
    Coordinate b_i = pinfo.getPos(newBond.i);
    Coordinate b_j = pinfo.getPos(newBond.j);
    
    std::vector <int> possibleBondNeighborTypes = {0,3};
    
    //sets initial bond neighbors of newly added bond
    pinfo.resetBondNeighborTagsByIndex(pinfo.getNumBonds()-1, possibleBondNeighborTypes);
    
    std::vector <int> newNeighboringTags = pinfo.getBondNeighborTagsByIndex(pinfo.getNumBonds()-1);
    
    // should measure distance between this new bond and the nucleator bead and if it is too close then remove it
    Coordinate nucleator_pos = pinfo.getPos(0);
    
    Coordinate dummy_k;
    
    
        
    Coordinate d_kl_nucleator = closestDistanceBetweenLineAndPoint(b_i, b_j, nucleator_pos, dummy_k);
    
    for(int neigh = newNeighboringTags.size()-1; neigh > -1; neigh--)
    {
        Bond neighborBond = pinfo.getBond(newNeighboringTags[neigh]);
        Coordinate neigh_i = pinfo.getPos(neighborBond.i);
        Coordinate neigh_j = pinfo.getPos(neighborBond.j);
        
        Coordinate dummy_l;

        Coordinate d_kl = closestDistanceBetweenLines(b_i, b_j, neigh_i, neigh_j, dummy_k, dummy_l);
        double currMinDist = d_kl.getMagnitude();
        if(currMinDist > rRepulsiveInteraction*5.0)
        {
            // remove this bond from newNeighboringTags
            newNeighboringTags[neigh] = newNeighboringTags.back();
            newNeighboringTags.pop_back();
        }
    }
    

    pinfo.setBondNeighborTagsByIndex(pinfo.getNumBonds()-1, newNeighboringTags);
    
    // append this bond to all of newNeighboringTags
    for(int i = 0; i < newNeighboringTags.size(); i++)
    {
        int neighboringBondTag = newNeighboringTags[i];
        
        pinfo.appendBondNeighborTag(neighboringBondTag, newBondTag);
    }
    
    return 0;
}




Coordinate Simulation::closestDistanceBetweenLineAndPoint(Coordinate a0, Coordinate a1, Coordinate b, Coordinate &aClose)
{
    //https://math.stackexchange.com/questions/1905533/find-perpendicular-distance-from-point-to-line-in-3d
    Coordinate a_disp = simBox.calcDisplacement(a1, a0);
    double a_disp_length = a_disp.getMagnitude();
    Coordinate a_unit = a_disp / a_disp_length;
    
    Coordinate ba_disp = simBox.calcDisplacement(b, a0);
    
    double t = ba_disp * a_unit;
    
    if(t < 0.0)
    {
        t = 0.0;
    }
    else if(t > a_disp_length)
    {
        t = a_disp_length;
    }
    
    aClose = a0 + t*a_unit;
    
    Coordinate d_kl = aClose - b;
    
    return d_kl;    
}


int Simulation::addFilament(int numBeads, Coordinate startPos, Coordinate direction, bool crosslinkNeighboring)
{   
    Filament newFil;
    
    
    newFil.setIsDaughter(false);
    filaments.push_back(newFil);
    
    int newFilamentTag = f_TaggedVector.add();
    
        
    int numStressMeasurements = 0;
    
    bool filamentBondedWithNeighbor = false;
    int numBondsWithNeighbors = 0;
    
    std::mt19937& threadEngine = generators[omp_get_thread_num()];
    
    std::vector <Bond> potentialBonds;
    
    for(int f = 0; f < numBeads; f++)
    {
        // should put new particle in grid and update neighbors of all associated particles
        Coordinate newPosition = {startPos.x + direction.x*L_zero*f, startPos.y + direction.y*L_zero*f, startPos.z + direction.z*L_zero*f};
        newPosition = simBox.periodicWrap(newPosition);
        Particle newParticle(newPosition.x, newPosition.y, newPosition.z);
        newParticle.setType(newFilamentTag);
        
        int newTag = addParticlePutInGrid(newParticle, numStressMeasurements);
        int priorTag = filaments[filaments.size()-1].addParticleBack(newTag);
        
        if(filamentBondsInBondList && priorTag > -1)
        {
            // adds a bond between priorTag and newTag of type 0 if they are on the same filament
            // as indicated from the type column in the xyz file
            Bond newBond = {priorTag, newTag, 0};
            newBond.creationTime = currSimTime;
            newBond.destructionTime = 1000000.0*140.0/40.0;
            
            int success_add = addBondFindNeighbors(newBond);
            
            if(success_add != 0)
            {
                // if this filament still exists, remove it
                if(f_TaggedVector.getIndexOfTag(newFilamentTag) >= 0)
                {
                    removeFilament(newFilamentTag);
                }
                
                return -1;
            }
        }
        
        // bond to existing filament beads that are within the correct range
        // bond if appropriate to existing particles
        
        
        if(crosslinkNeighboring == true && numBondsWithNeighbors < 5)
        {
            Coordinate particleJ = pinfo.getPosByIndex(pinfo.getNumParticles()-1);
            int tempBondJTag = newTag;
        
            // this should use the grid rather than an exhaustive search
            for(int loopI = 0; loopI < pinfo.getNumParticles()-1-f; loopI++)
            {
                Coordinate particleI = pinfo.getPosByIndex(loopI);
                int tempBondITag = pinfo.getTagAtIndex(loopI);

                Coordinate uij = particleI - particleJ;
                
                uij = simBox.periodicWrap(uij);
            
                double distMag = uij.getMagnitude();
                
                if(distMag >= minCrosslinkBondDist && distMag <= maxCrosslinkBondDist)
                {
                    {
                        // assumes that crosslinker bonds are stored as type 1 in bondTypes
                        Bond newBond = {tempBondITag, tempBondJTag, 1};
                        newBond.creationTime = currSimTime;
                        newBond.destructionTime = 1000000.0*140.0/40.0;
                        
                        bool bondAlreadyExists = pinfo.checkBondExists(newBond);
                        
                        if(bondAlreadyExists == false)
                        {
                            potentialBonds.push_back(newBond);
                        }
                    }
                }
            }
            
            
        }
    }
    
    return newFilamentTag;
        
}


int Simulation::addFilamentSpecifyBondType(int numBeads, Coordinate startPos, Coordinate direction, bool crosslinkNeighboring, int bondIndex)
{   
    Filament newFil;
    newFil.setIsDaughter(false);
    filaments.push_back(newFil);
    int newFilamentTag = f_TaggedVector.add();
        
    int numStressMeasurements = 0;
    
    bool filamentBondedWithNeighbor = false;
    int numBondsWithNeighbors = 0;
    
    std::mt19937& threadEngine = generators[omp_get_thread_num()];
    
    std::vector <Bond> potentialBonds;
    
    double bondLength = bondTypes[bondIndex].getEqDist();
    
    
    for(int f = 0; f < numBeads; f++)
    {
        // should put new particle in grid and update neighbors of all associated particles
        Coordinate newPosition = {startPos.x + direction.x*bondLength*f, startPos.y + direction.y*bondLength*f, startPos.z + direction.z*bondLength*f};
        newPosition = simBox.periodicWrap(newPosition);
        Particle newParticle(newPosition.x, newPosition.y, newPosition.z);
        newParticle.setType(newFilamentTag);
        
        int newTag = addParticlePutInGrid(newParticle, numStressMeasurements);
        int priorTag = filaments[filaments.size()-1].addParticleBack(newTag);
        
        if(filamentBondsInBondList && priorTag > -1)
        {
            // adds a bond between priorTag and newTag of type 0 if they are on the same filament
            // as indicated from the type column in the xyz file
            Bond newBond = {priorTag, newTag, bondIndex};
            newBond.creationTime = currSimTime;
            newBond.destructionTime = 1000000.0*140.0/40.0;
            
            int success_add = addBondFindNeighbors(newBond);
            
            if(success_add != 0)
            {
                // if this filament still exists, remove it
                if(f_TaggedVector.getIndexOfTag(newFilamentTag) >= 0)
                {
                    removeFilament(newFilamentTag);
                }
                
                return -1;
            }
        }
        
        // bond to existing filament beads that are within the correct range
        // bond if appropriate to existing particles
        
        
        if(crosslinkNeighboring == true && numBondsWithNeighbors < 5)
        {
            Coordinate particleJ = pinfo.getPosByIndex(pinfo.getNumParticles()-1);
            int tempBondJTag = newTag;
        
            // this should use the grid rather than an exhaustive search
            for(int loopI = 0; loopI < pinfo.getNumParticles()-1-f; loopI++)
            {
                Coordinate particleI = pinfo.getPosByIndex(loopI);
                int tempBondITag = pinfo.getTagAtIndex(loopI);

                Coordinate uij = particleI - particleJ;
                
                uij = simBox.periodicWrap(uij);
            
                double distMag = uij.getMagnitude();
                
                if(distMag >= minCrosslinkBondDist && distMag <= maxCrosslinkBondDist)
                {
                    // peramanent bond
                    {
                        // assumes that crosslinker bonds are stored as type 1 in bondTypes
                        Bond newBond = {tempBondITag, tempBondJTag, 1};
                        newBond.creationTime = currSimTime;
                        // never destruct permanent bonds
                        newBond.destructionTime = currSimTime + 10.0*140.0/40.0;
                        
                        bool bondAlreadyExists = pinfo.checkBondExists(newBond);
                        
                        if(bondAlreadyExists == false)
                        {
                            potentialBonds.push_back(newBond);
                        }
                    }
                }
            }
            
            
        }
    }
    
    
    for(int i = 0; i < potentialBonds.size(); i++)
    {
        if(i >= 5)
        {
            // max 5 bonds
            break;
        }
        
        int rndIndex = (int)(potentialBonds.size()*dis(threadEngine));
        Bond rndBond = potentialBonds[rndIndex];
        
        pinfo.addBond(rndBond);
        
        potentialBonds[rndIndex] = potentialBonds.back();
        potentialBonds.pop_back();
    }
    
    return newFilamentTag;
}

void Simulation::addRandomFilament(bool crosslinkNeighboring)
{
    double randX, randY, randZ;
    
    std::mt19937& threadEngine = generators[omp_get_thread_num()];
    
    double rndVal = dis(threadEngine);
    
    randZ = -1.0;
    randY = dis(threadEngine)*simBox.getBoxLength(1) - 0.5*simBox.getBoxLength(1);
    
    double dirX, dirY, dirZ;
    
    if(rndVal < 0.50)
    {
        double theta = pi*0.5;
        
        randX = dis(threadEngine)*simBox.getBoxLength(0) - 0.5*simBox.getBoxLength(0);
        
        double randomAngle = theta + 2.44346*dis(threadEngine) - 1.22173;
        double randomOutOfPlaneAngle = 10.0*pi/180.0 * 2.0 * (dis(threadEngine) - 0.5);
        
        dirX = cos(randomOutOfPlaneAngle)*cos(randomAngle);
        dirZ = cos(randomOutOfPlaneAngle)*sin(randomAngle);
        dirY = sin(randomOutOfPlaneAngle);
    }
    else
    {
        randX = dis(threadEngine)*simBox.getBoxLength(0) - 0.5*simBox.getBoxLength(0);
        
        double randomOutOfPlaneAngle = 10.0*pi/180.0 * 2.0 * (dis(threadEngine) - 0.5);
        
        dirX = 0.0;
        dirZ = 1.0;
        dirY = 0.0;
    }
    
    if(dirZ < 0.0)
    {
        dirX = copysign(1.0, dirX);
        dirZ = 0.0;
    }
    
    int filamentLength = (int)round(1.0 / L_zero) + 1;
    
    Coordinate firstBead = {randX, randY, randZ};
    Coordinate direction = {dirX, dirY, dirZ};
    
    addFilament(filamentLength, firstBead, direction, crosslinkNeighboring);
}


void Simulation::addFilamentTouchingNucleatorBead(bool crosslinkNeighboring)
{
    // assumes that the nucleator bead is always the first bead in pinfo
    Coordinate nucleatorBeadPos = pinfo.getPosByIndex(0);
    
    std::mt19937& threadEngine = generators[omp_get_thread_num()];
    
    double randU = 2.0*dis(threadEngine) - 1.0;
    double randTheta = 2.0*pi*dis(threadEngine);
    
    Coordinate randUnitVector = Coordinate {sqrt(1.0 - randU*randU) * cos(randTheta), sqrt(1.0 - randU*randU) * sin(randTheta), randU};
    
    // this is barbed end position
    Coordinate surfaceTouchPoint = nucleatorBeadPos + nucleatorBeadRadius * randUnitVector;
    
    double randU_dir = 2.0*dis(threadEngine) - 1.0;
    double randTheta_dir = 2.0*pi*dis(threadEngine);
    
    Coordinate randUnitVector_dir = Coordinate {sqrt(1.0 - randU_dir*randU_dir) * cos(randTheta_dir), sqrt(1.0 - randU_dir*randU_dir) * sin(randTheta_dir), randU_dir};
    
    Coordinate direction = randUnitVector_dir;
    
    int filamentNumBeads = 2;
    
    Coordinate pointedEnd;
    
    if(direction*randUnitVector > 0.0)
    {
        //then direction is pointing outside the sphere
        pointedEnd = surfaceTouchPoint;
    }
    else
    {
        //then direction is pointing into the sphere
        pointedEnd = surfaceTouchPoint + 2.7E-3*(starting_bond_type-3) * (-direction) * (filamentNumBeads - 1);
    }
    
    addFilamentSpecifyBondType(filamentNumBeads, pointedEnd, direction, crosslinkNeighboring, starting_bond_type);
}

// branch off of any bead in tagsPotentialCrosslink which does not already have a branch
void Simulation::addRandomBranch(int numBeads)
{
    std::vector<int> tagsPotentialCrosslink = particleGrid.getTagsInRegion(pinfo.getPosByIndex(0), nucleatorBeadRadius*1.25);
    
    for(int i = tagsPotentialCrosslink.size()-1; i > -1 ; i--)
    {
        double dist = simBox.calcDistance(pinfo.getPos(tagsPotentialCrosslink[i]), pinfo.getPosByIndex(0));
        
        if(((dist - nucleatorBeadRadius) > 0.1) || tagsPotentialCrosslink[i] == 0)
        {
            // remove this tag from tagsPotentialCrosslink
            tagsPotentialCrosslink[i] = tagsPotentialCrosslink.back();
            tagsPotentialCrosslink.pop_back();
        }
    }

    std::mt19937& threadEngine = generators[omp_get_thread_num()];  

    bool notYetAdded = true;
    
    //#pragma omp parallel for schedule(static, 1)
    for(int curr_pos = 0; curr_pos < tagsPotentialCrosslink.size(); curr_pos++)
    
    {
       
        int randomTag = tagsPotentialCrosslink[curr_pos];
        
        bool skipBead = false;
        // check that there is not a crosslink already off of randomTag
        std::vector <struct HalfBond> tempBondList = pinfo.getBondTags(randomTag);
        
        
        for(int bb = 0; bb < tempBondList.size(); bb++)
        {
            Bond tempBond = pinfo.getBond(tempBondList[bb].bondTag);
            
            if(tempBond.bondTypeIndex == 3)
            {
                // there is already a branch off of this particle, so do not bond to it
                skipBead = true;
                break;
            }
        }
        
        Coordinate beadPos = pinfo.getPos(randomTag);
        
        double dist = simBox.calcDistance(beadPos, pinfo.getPosByIndex(0));
        
        int filamentTag = pinfo.getBeadToFilament(randomTag);
        int filamentIndex = f_TaggedVector.getIndexOfTag(filamentTag);
        
        int filamentPosRandomBead = -1;
        // need to use filaments vector since branch beads arent guaranteed to be sequential
        for(int f = 0; f < filaments[filamentIndex].getNumParticles(); f++)
        {
            if(filaments[filamentIndex].getTagAtIndex(f) == randomTag)
            {
                filamentPosRandomBead = f;
                break;
            }
        }
        
        int bondedTag = -1;
        
        // filamentPosRandomBead can't be the barbed end since there will not be a third bead for branch angle
        if(filamentPosRandomBead < filaments[filamentIndex].getNumParticles()-1)
        {
            bondedTag = filaments[filamentIndex].getTagAtIndex(filamentPosRandomBead + 1);
        }
        else
        {
            skipBead = true;
        }
        if(!skipBead && pinfo.checkTagExists(bondedTag))
        {
            int bondedFilamentTag = pinfo.getBeadToFilament(bondedTag);
           
            // if this is on the same filament then try to add crosslink

            double randomRoll = dis(threadEngine);
            double Rdt = 1.0 - exp(-branchingRatePerSegment*timeBetweenBranchingEvents);
            
            if(bondedFilamentTag == filamentTag && randomRoll < Rdt)
            {
                // 70 degrees in radians
                double theta = 1.22173;
                
                double u = cos(theta);
                
                
                double phi = dis(threadEngine) * 2.0 * pi;
                
                
                Coordinate randPositionOnConeUnit = Coordinate {sqrt(1.0-u*u)*cos(phi), sqrt(1.0-u*u)*sin(phi), u};
                

                Coordinate bondedPos = pinfo.getPos(bondedTag);
                
                Coordinate filamentUnitVector = (bondedPos - beadPos).getUnitCoord();
                
                
                // now need to rotate North pole of cone to angle that actin filament is at
                double initialTheta = acos(filamentUnitVector.z / 1.0);
                double initialPhi = atan2(filamentUnitVector.y, filamentUnitVector.x);
                

                struct Coordinate rotK = Coordinate {0.0, 1.0, 0.0};
                
                // rotate by initialTheta
                
                double currDotProduct = rotK*randPositionOnConeUnit;
                struct Coordinate currCrossProduct = rotK.crossProduct(randPositionOnConeUnit);
                
                randPositionOnConeUnit = randPositionOnConeUnit*cos(initialTheta) + currCrossProduct*sin(initialTheta) + rotK*currDotProduct*(1.0 - cos(initialTheta));
                
                rotK.x = 0.0;
                rotK.y = 0.0;
                rotK.z = 1.0;
                
                // rotate by initialPhi
                currDotProduct = rotK*randPositionOnConeUnit;
                currCrossProduct = rotK.crossProduct(randPositionOnConeUnit);
                
                randPositionOnConeUnit = randPositionOnConeUnit*cos(initialPhi) + currCrossProduct*sin(initialPhi) + rotK*currDotProduct*(1.0 - cos(initialPhi));
    
                Coordinate direction = randPositionOnConeUnit;

                
                // keep this at currNumParticles so that the pointed tag is at pointed end
                int currNumParticles = pinfo.getNumParticles();

                // number of particles in filament is first argument
                
                double initialLength = bondTypes[3].getEqDist();
                int new_filament_tag = addFilamentSpecifyBondType(numBeads, beadPos + initialLength * direction, direction, false, starting_bond_type);
                
                if(new_filament_tag != -1)
                {
                    // bond between pointed end and randomTag
                    Bond newBond = {randomTag, pinfo.getTagAtIndex(currNumParticles), 3};
                    
                    newBond.creationTime = currSimTime;
                    newBond.destructionTime = 100000.0*140.0/40.0;

                    // if this returns anythign other than 0 need to remove all bonds for this branch
                    int success_add = addBondFindNeighbors(newBond);
                    

                    if(success_add != 0)
                    {
                        removeFilament(new_filament_tag);
                    }
                    else
                    {
                        // angle between barbed end, randomTag and bondedTag (70 degrees)
                        Angle newAngle = {bondedTag, randomTag, pinfo.getTagAtIndex(currNumParticles), 0};
                        pinfo.addAngle(newAngle);
                        
                        // straight angle
                        newAngle = {randomTag, pinfo.getTagAtIndex(currNumParticles), pinfo.getTagAtIndex(currNumParticles+1), 1};
                        pinfo.addAngle(newAngle);
                    }
                }
            }
        }
        
    }    
}

void Simulation::removeFilament(int filamentTag)
{
    int filamentIndex = f_TaggedVector.getIndexOfTag(filamentTag);
    
    int initialNumParticles = filaments[filamentIndex].getNumParticles();
    for(int i = 0; i < initialNumParticles; i++)
    {        
        int particleTag = filaments[filamentIndex].getBarbedTag();
        filaments[filamentIndex].removeParticleBack();
        
        pinfo.removeParticle(particleTag);
        particleGrid.removeFromGrid(particleTag);
    }
    
    filaments[filamentIndex] = filaments.back();
    filaments.pop_back();
    f_TaggedVector.remove(filamentTag);
}

void Simulation::removeFilamentsNearPos(Coordinate pos, double dist)
{
    int initialNumParticles = pinfo.getNumParticles();
    std::vector <int> removeFilamentsVector;
    
    // if a particle on a filament is near the focal adhesion mark it for removal via removeFilament function
    // should only search through particles which are near the focal adhesion
    #pragma omp parallel for
    for(int p = 0; p < pinfo.getNumParticles(); p++)
    {
        Coordinate currPos = pinfo.getPosByIndex(p);
        
        if(currPos.z <= pos.z + dist && currPos.z >= pos.z - dist && currPos.y <= 0.0)
        {
            Coordinate currDisp = simBox.calcDisplacement(currPos, pos);
            currPos.y = 0.0;
            
            double currDist = currDisp.getMagnitude();
            if(currDist <= dist)
            {
                int currFilamentTag = pinfo.getBeadToFilamentByIndex(p);
                
                if(currFilamentTag >= 0)
                {
                    #pragma omp critical
                    {
                        bool found = false;
                        
                        for(int m = 0; m < removeFilamentsVector.size(); m++)
                        {
                            if(removeFilamentsVector[m] == currFilamentTag)
                                found = true;
                        }
                        
                        if(found == false)
                            removeFilamentsVector.push_back(currFilamentTag);
                    }
                }
            }
        }
    }
    
    for(int f = 0; f < removeFilamentsVector.size(); f++)
    {
        removeFilament(removeFilamentsVector[f]);
    }
}


void Simulation::removeAllCrosslinks()
{
    for(int b = pinfo.getNumBonds()-1; b > -1; b--)
    {
        Bond tempBond = pinfo.getBondByIndex(b);
        
        // crosslinks have a type of 1 or 2
        if(tempBond.bondTypeIndex != 0)
        {
            pinfo.removeBondByIndex(b);
        }
    }
}

void Simulation::removeTemporaryCrosslinks()
{
    for(int b = pinfo.getNumBonds()-1; b > -1; b--)
    {
        Bond tempBond = pinfo.getBondByIndex(b);
        
        // temporary crosslinks are type 2
        if(tempBond.bondTypeIndex == 2)
        {
            pinfo.removeBondByIndex(b);
        }
    }
}

void Simulation::removeBond(int bondTag)
{
    int bondIndex = pinfo.getBondPosOfTag(bondTag);
    this->removeBondByIndex(bondIndex);
}

void Simulation::removeBondByIndex(int bondIndex)
{
    Bond currBond = pinfo.getBondByIndex(bondIndex);
    
    int particleITag = currBond.i;
    int particleJTag = currBond.j;
    
    int filamentI = pinfo.getBeadToFilament(particleITag);
    int filamentJ = pinfo.getBeadToFilament(particleJTag);
    
    pinfo.removeBondByIndex(bondIndex);

    if(filamentI == filamentJ)
    {
        // this bond is internal to a filament, need to break this filament at this bond
        
        int filamentIndexI = f_TaggedVector.getIndexOfTag(filamentI);
        
        int indexI = filaments[filamentIndexI].getIndexOfTag(particleITag);
        int indexJ = filaments[filamentIndexI].getIndexOfTag(particleJTag);
        
        if(indexI > indexJ)
        {
            int tempSwap = indexI;
            indexI = indexJ;
            indexJ = tempSwap;
        }
        
        Filament tempFilament = filaments[filamentIndexI].breakFilamentAtPos(indexI);
        
        cleanExistingFilament(filamentIndexI, filamentI);
        cleanNewFilament(tempFilament);
    }
}

void Simulation::removeParticlesOutsideBox()
{
    double minZ = simBox.getBoxLength(2);

    for(int p = pinfo.getNumParticles()-1; p > -1; p--)
    {
        int tempTag = pinfo.getTagAtIndex(p);
        
        if(tempTag >= 0 && pinfo.getPosByIndex(p).z < -minZ)
        {
            removeParticleByIndex(p);
        }
    }
}

void Simulation::cleanExistingFilament(int filamentIndex, int filamentTag)
{
    if(filaments[filamentIndex].getNumParticles() == 1)
    {
        // only one particle in filament, remove it and the filament
        int onlyTag = filaments[filamentIndex].getTagAtIndex(0);
        
        // removes any bonds associated with this particle as well in ParticleInfo.cpp
        pinfo.removeParticle(onlyTag);
        particleGrid.removeFromGrid(onlyTag);

        filaments[filamentIndex] = filaments.back();
        filaments.pop_back();
        f_TaggedVector.remove(filamentTag);
    }
    else if(filaments[filamentIndex].getNumParticles() == 0)
    {
        // no particles in filament, remove filament
        filaments[filamentIndex] = filaments.back();
        filaments.pop_back();
        f_TaggedVector.remove(filamentTag);
    }
}

void Simulation::cleanNewFilament(Filament tempFil)
{
    if(tempFil.getNumParticles() >= 2)
    {
        filaments.push_back(tempFil);
        
        // don't cap newly broken filaments
        filaments[filaments.size()-1].setCappedBarbedEnd(false);
        int newFilamentTag = f_TaggedVector.add();
        
        // update which filament beads are in for the new filament
        for(int v = 0; v < tempFil.getNumParticles(); v++)
        {
            int currTag = tempFil.getTagAtIndex(v);
            
            pinfo.setBeadToFilament(currTag, newFilamentTag);
        }
    }
    else if(tempFil.getNumParticles() == 1)
    {
        int onlyTag = tempFil.getTagAtIndex(0);
        
        // always remove regardless of any bonds to onlyTag
        pinfo.removeParticle(onlyTag);
        particleGrid.removeFromGrid(onlyTag);
    }
}


void Simulation::removeParticleByIndex(int pIndex)
{
    int pTag = pinfo.getTagAtIndex(pIndex);
    
    // need to update filament indexing
    int filamentTag = pinfo.getBeadToFilamentByIndex(pIndex);
    
    Coordinate tempPos = pinfo.getPos(pTag);

    if(filamentTag >= 0)
    {
        int filamentIndex = f_TaggedVector.getIndexOfTag(filamentTag);
        
        int filamentSize = filaments[filamentIndex].getNumParticles();
        
        if(filamentSize > 2)
        {        
            Filament tempFil = filaments[filamentIndex].removeBeadByTag(pTag);
            cleanExistingFilament(filamentIndex, filamentTag);
            // adds tempFil to filaments if it is big enough
            cleanNewFilament(tempFil);
        }
        else if(filamentSize == 2)
        {
            // filament with bead to be removed it only size 2, remove both bead and filament
            
            int otherTag = filaments[filamentIndex].getTagAtIndex(0);
            
            if(otherTag == pTag)
            {
                otherTag = filaments[filamentIndex].getTagAtIndex(1);
            }
            
            filaments[filamentIndex] = filaments.back();
            filaments.pop_back();
            f_TaggedVector.remove(filamentTag);
            
            // if this tag is of a lower index it may have not been given a bin in the grid yet by putAllParticlesInGrid
            pinfo.removeParticle(otherTag);
            particleGrid.removeFromGrid(otherTag);
        }
        else
        {
            // filament with bead to be removed is size 1 (or 0 if there is an error), remove filament
            
            filaments[filamentIndex] = filaments.back();
            filaments.pop_back();
            f_TaggedVector.remove(filamentTag);
        }
        
        
        // actually remove the particle from pinfo (also removes the bonds (but not yet angles))
        pinfo.removeParticle(pTag);
        particleGrid.removeFromGrid(pTag);
    }
}

void Simulation::addBeadCenterBox()
{
     // nucleator bead
    Particle nucleatorBead(0.0, 0.0, 0.0);
    nucleatorBead.setType(-1);
    int numStressMeasurements = 0;
    int newTag = pinfo.addParticle(nucleatorBead, numStressMeasurements);
}

void Simulation::addSegmentToExistingFilament(int filamentIndex)
{
    std::mt19937& threadEngine = generators[omp_get_thread_num()];  

    if(filaments[filamentIndex].getCappedBarbedEnd() == false)
    {
        Coordinate barbedEndPosition = pinfo.getPos(filaments[filamentIndex].getBarbedTag());
         
        Coordinate nextPosition = pinfo.getPos(filaments[filamentIndex].getTagAtIndex(filaments[filamentIndex].getNumParticles()-2));       
        
        Coordinate unitVector = simBox.calcDisplacement(barbedEndPosition, nextPosition).getUnitCoord();
        
        Coordinate newPosition = barbedEndPosition + unitVector*L_zero;
        
        Coordinate nucleatorBeadPos = pinfo.getPosByIndex(0);
        
        double tmpDist = simBox.calcDistance(newPosition, nucleatorBeadPos);
        
        if(filaments[filamentIndex].getNumParticles() >= 3)
        {
            // cap this filament
            {
                filaments[filamentIndex].setCappedBarbedEnd(true);
            }
        }
        
        if(filaments[filamentIndex].getCappedBarbedEnd() == false)
        {
            int filamentTag = f_TaggedVector.getTagAtIndex(filamentIndex);
            
            newPosition = simBox.periodicWrap(newPosition);
            Particle newParticle(newPosition.x, newPosition.y, newPosition.z);
            newParticle.setType(filamentTag);

            #pragma omp critical
            {
                int newTag = addParticlePutInGrid(newParticle, 0);
                int priorTag = filaments[filamentIndex].addParticleBack(newTag);
                
                if(filamentBondsInBondList && priorTag > -1)
                {
                    // adds a bond between priorTag and newTag of type 0 if they are on the same filament
                    // as indicated from the type column in the xyz file
                    Bond newBond = {priorTag, newTag, 0};
                    newBond.creationTime = currSimTime;
                    newBond.destructionTime = 1000000.0*140.0/40.0;
                    
                    // okay to not check return value since only adding a single bond
                    addBondFindNeighbors(newBond);
                }
                else
                {
                    cout << "error filament: " << filamentIndex << " " << filamentTag << endl;
                    exit(0);
                }
            }
        }
    }
}


void Simulation::addMonomerToExistingFilament(int filamentIndex)
{
    std::mt19937& threadEngine = generators[omp_get_thread_num()];      
    int numParticlesInFilament = filaments[filamentIndex].getNumParticles();
    
    if(filaments[filamentIndex].getCappedBarbedEnd() == false)
    {
        if(filaments[filamentIndex].getCappedBarbedEnd() == false)
        {
            int barbedEndTag = filaments[filamentIndex].getTagAtIndex(numParticlesInFilament-1);
            int nextTag = filaments[filamentIndex].getTagAtIndex(numParticlesInFilament-2);
            
            
            int bondTag = pinfo.findBondTag(nextTag, barbedEndTag);
            
            int bondType = pinfo.getBond(bondTag).bondTypeIndex;
            
            Coordinate barbedEndPosition = pinfo.getPos(barbedEndTag);
            Coordinate nextPosition = pinfo.getPos(nextTag);
            
            Coordinate unit_disp = (simBox.calcDisplacement(barbedEndPosition, nextPosition)).getUnitCoord();
            
            if(bondType >= 4)
            {
                int newBondType = bondType + 1;
                
                // also need to move particle at free end here to simulate growth rather than expansion of filament bond both directions
                
                pinfo.setPos(barbedEndTag, nextPosition + bondTypes[newBondType].getEqDist()*unit_disp);
                
                if(newBondType == 2+numMonomersPerBead)
                {
                    // switch this to type 0
                    newBondType = 0;
                    
                    
                    if(numParticlesInFilament >= maxLength)
                    {
                        // cap this filament
                        
                        {
                            filaments[filamentIndex].setCappedBarbedEnd(true);
                        }
                    }
                }
                
                pinfo.setBondType(bondTag, newBondType);
            }
            else if(bondType == 0)
            {
                
                if(numParticlesInFilament >= maxLength)
                {
                    filaments[filamentIndex].setCappedBarbedEnd(true);
                }
                else
                {
                    Coordinate newPosition = barbedEndPosition + 2.7E-3*(starting_bond_type-3)*unit_disp;
                    
                    Coordinate nucleatorBeadPos = pinfo.getPosByIndex(0);
            
            
                    int filamentTag = f_TaggedVector.getTagAtIndex(filamentIndex);
                    // add a new bead then with a bond of type 4

                    newPosition = simBox.periodicWrap(newPosition);
                    Particle newParticle(newPosition.x, newPosition.y, newPosition.z);
                    newParticle.setType(filamentTag);

                    #pragma omp critical
                    {
                        int newTag = addParticlePutInGrid(newParticle, 0);
                        int priorTag = filaments[filamentIndex].addParticleBack(newTag);
                        

                        if(filamentBondsInBondList && priorTag > -1)
                        {
                            // adds a bond between priorTag and newTag of type 0 if they are on the same filament
                            // as indicated from the type column in the xyz file
                            Bond newBond = {priorTag, newTag, starting_bond_type};
                            newBond.creationTime = currSimTime;
                            
                            newBond.destructionTime = 1000000.0*140.0/40.0;
                            
                            addBondFindNeighborsNoRemoval(newBond);
                        }
                        else
                        {
                            cout << "error filament: " << filamentIndex << " " << filamentTag << endl;
                            exit(0);
                        }
                    }
                }
            }
        }
    }
}

void Simulation::addSegmentsToAllFilaments(double rate, double time_between_attempts)
{
    //#pragma omp parallel for
    for(int f = 0; f < filaments.size(); f++)
    {
        std::mt19937& threadEngine = generators[omp_get_thread_num()]; 
        
        double randVal = dis(threadEngine);
        
        
        // measure force on bond closest to barbed end of this filament
        // filament should always have at least 2 beads so this should be safe
        
        if(filaments[f].getNumParticles() < 2)
        {
            std::cout << "Not enough particles in " << f << " " << filaments[f].getNumParticles() << std::endl;
        }
        
        int barbedEndTag = filaments[f].getBarbedTag();
        int neighboringBarbedEndTag = filaments[f].getTagAtIndex(filaments[f].getNumParticles()-2);
        
        Coordinate barbedEndTagPos = pinfo.getPos(barbedEndTag);
        Coordinate neighboringTagPos = pinfo.getPos(neighboringBarbedEndTag);
        
        Coordinate bc = simBox.calcDisplacement(barbedEndTagPos, neighboringTagPos);
        double bcMag = bc.getMagnitude();
        
        int tempBondTag = pinfo.findBondTag(barbedEndTag, neighboringBarbedEndTag);
        Bond tempBond = pinfo.getBond(tempBondTag);
        double bondForce = bondTypes[tempBond.bondTypeIndex].calcForce(bcMag);
        
        double rate_multiplier = 1.0;
        
        
        if(bondForce > 0.0)
        {
            // bond is under compression, Brownian ratchet expression from Mogilner and Oster 1996
            rate_multiplier = exp(-bondForce * 2.7 / (4.114));
        }
        
        double Rdt = rate * rate_multiplier * time_between_attempts;
        
        if(randVal < 1.0 - exp(-Rdt))
        {
            //#pragma omp critical
            {
                addMonomerToExistingFilament(f);
            }
        }
    }
}

Coordinate Simulation::closestDistanceBetweenLines(Coordinate a0, Coordinate a1, Coordinate b0, Coordinate b1, Coordinate &aClose, Coordinate &bClose)
{
    //Based on https://stackoverflow.com/questions/2824478/shortest-distance-between-two-line-segments

    // Given two lines defined by (a0,a1,b0,b1)
    // Return the closest points on each segment and their distance
    //

    Coordinate d_kl;

    // Calculate denomitator
    Coordinate A = a1 - a0;
    Coordinate B = b1 - b0;
    double magA = A.getMagnitude();
    double magB = B.getMagnitude();

    Coordinate _A = A / magA;
    Coordinate _B = B / magB;

    Coordinate cross = _A.crossProduct(_B);
    double denom = cross.getMagnitude() * cross.getMagnitude();


    // If lines are parallel (denom=0) test if lines overlap.
    // If they don't overlap then there is a closest point solution.
    // If they do overlap, there are infinite closest positions, but there is a closest distance
    if (denom == 0)
    {
        double d0 = _A * (b0-a0);

        // Overlap only possible with clamping

        double d1 = _A * (b1-a0);

        if (d0 <= 0 && 0 >= d1)
        {
            if (fabs(d0) < fabs(d1))
            {
                aClose = a0;
                bClose = b0;
                
                d_kl = (a0 - b0);
                
                //tk = 0.5;
                //tl = 0.5;

                return d_kl;
            }
            aClose = a0;
            bClose = b1;
            
            //tk = 0.5;
            //tl = -0.5;
            d_kl = (a0 - b1);

            return d_kl;
        }
        else if (d0 >= magA && magA <= d1)
        {
            if (fabs(d0) < fabs(d1))
            {
                // tl is -0.5, tk is 0.5
                aClose = a1;
                bClose = b0;

                d_kl = (a1 - b0);
                
                return d_kl;
            }
            aClose = a1;
            bClose = b1;
                
            d_kl = (a1 - b1);
                
            return d_kl;
        }

        // Segments overlap, return distance between parallel segments
        aClose = a0 + A*0.5;
        bClose = b0 + B*0.5;
            
        d_kl = (((d0 * _A) + a0) - b0);
        return d_kl;
    }


    // Lines criss-cross: Calculate the projected closest points
    Coordinate t = (b0 - a0);
    double detA = calcDeterminant(t, _B, cross);
    double detB = calcDeterminant(t, _A, cross);

    double t0 = detA / denom;
    double t1 = detB / denom;

    Coordinate pA = a0 + (_A * t0); // Projected closest point on segment A
    Coordinate pB = b0 + (_B * t1); // Projected closest point on segment B


    // Clamp projections
    if (t0 < 0)
        pA = a0;
    else if (t0 > magA)
        pA = a1;

    if (t1 < 0)
        pB = b0;
    else if (t1 > magB)
        pB = b1;

    double dot;
    // Clamp projection A
    if (t0 < 0 || t0 > magA)
    {
        dot = _B * (pA - b0);
        if (dot < 0)
            dot = 0;
        else if (dot > magB)
            dot = magB;
        pB = b0 + (_B * dot);
    }
    // Clamp projection B
    if (t1 < 0 || t1 > magB)
    {
        dot = _A*(pB - a0);
        if (dot < 0)
            dot = 0;
        else if (dot > magA)
            dot = magA;
        pA = a0 + (_A * dot);
    }
    aClose = pA;
    bClose = pB;
    
    d_kl = (pA - pB);
    return d_kl;
}

double Simulation::calcDeterminant(Coordinate firstRow, Coordinate secondRow, Coordinate thirdRow)
{
    double x = secondRow.y*thirdRow.z - thirdRow.y*secondRow.z;
    double y = secondRow.x*thirdRow.z - thirdRow.x*secondRow.z;
    double z = secondRow.x*thirdRow.y - thirdRow.x*secondRow.y;
    
    return firstRow.x*x - firstRow.y*y + firstRow.z*z;
}



void Simulation::run()
{
    bool verbose = false;
    
    ofstream outputXYZ;
    outputXYZ.open ("out-" + outputFileName + ".xyz");
    
    ofstream outputBND;
    outputBND.open ("out-" + outputFileName + ".bnd");
    
    ofstream outputANG;
    outputANG.open ("out-" + outputFileName + ".ang");
    
    double t1,t2;
	t1=omp_get_wtime();
    
    double initialSimTime = currSimTime;

    omp_set_num_threads(numThreads);
    
    for(int i = 0; i < omp_get_max_threads(); i++)
    {
        vector <Coordinate> tempVec;
        
        for(int j = 0; j < pinfo.getNumParticles(); j++)
        {
            Coordinate c = {0.0, 0.0, 0.0};
            tempVec.push_back(c);
        }
        tempForces.push_back(tempVec);
    }
    
    double timeBetweenAddFilaments = 0.01*140.0/40.0;
    int timestepAddFilaments = (int)round(timeBetweenAddFilaments / dt);
    cout << "timestepAddFilaments: " << timestepAddFilaments << endl;
    int addFilamentCount = 0;
    
    timeBetweenBranchingEvents = 0.001*(140.0/40.0);
    int timestepAddBranch = (int)round(timeBetweenBranchingEvents / dt);
    timeBetweenBranchingEvents = timestepAddBranch*dt;
    
    int numMonomersPerSecondAdd = 40;
    double numSegmentsPerSecondAdd = (double)numMonomersPerSecondAdd / numMonomersPerBead;
    double timeBetweenGrowingFilaments = dt;
    int timestepGrowFilaments = (int)round(timeBetweenGrowingFilaments / dt);
    timeBetweenGrowingFilaments = timestepGrowFilaments * dt;
    
    double scaled_denovo_rate_real = 1.131768484 * 4*pi*nucleatorBeadRadius*nucleatorBeadRadius;
    
    double time_between_denovo = 1.0 / scaled_denovo_rate_real * 40.0 / 140.0 * (140.0/40.0);
    int timestep_denovo = (int)round(time_between_denovo / dt);
    time_between_denovo = timestep_denovo * dt;
    
    
    
    double lastSnapshotTime = omp_get_wtime();
    double currSnapshotTime = omp_get_wtime();
    
    if(pinfo.getNumParticles() == 0)
    {
        addBeadCenterBox();
    }
    
    int count = 1;
    
    for(int mainLoopVar = 0; mainLoopVar < numSteps; mainLoopVar++)
    {
        if(verbose == true)
            cout << "before putting in grid" << endl;
        
        
        if(mainLoopVar % updateGridStep == 0)
        {
            // remove particles here if they are outside the simulation box when check if they are in grid
            // also generate the vector for use with addCrosslinks function
            putAllParticlesInGrid();
        }
        
        // calculate main forces
        double leadingEdgeForce = calcForces();
        
        
        // take a snapshot
        if(mainLoopVar % snapshotStep == 0)
        {            
            currSnapshotTime = omp_get_wtime();
            
            
            cout << fixed << currSimTime << " / " << (totalSimulationTime + initialSimTime) << " time(s): " << currSnapshotTime - lastSnapshotTime << endl;
            
            lastSnapshotTime = currSnapshotTime;
            
            writeXYZFile(outputXYZ);
            writeBondFile(outputBND);
            writeAngleFile(outputANG);
            
            ofstream outputLastXYZ;
            outputLastXYZ.open ("last-" + outputFileName + ".xyz");
            
            ofstream outputLastBND;
            outputLastBND.open ("last-" + outputFileName + ".bnd");
            
            ofstream outputLastANG;
            outputLastANG.open ("last-" + outputFileName + ".ang");
            
            writeXYZFile(outputLastXYZ);
            writeBondFile(outputLastBND);
            writeAngleFile(outputLastANG);
            
            outputLastXYZ.close();
            outputLastBND.close();
            outputLastANG.close();
        }       
        
        moveParticles();
        
        currSimTime += dt;
        
    
        // update actin network connectivity here
        
        int numRemovedBonds = removeFlaggedBonds();
        
      
        
        if(pinfo.getNumParticles() == 1)
        {
            for(int i = 0; i < 200; i++)
            {
                addFilamentTouchingNucleatorBead(false);
            }
        }
        
        
        // add a random filament to the network if it is time to do so
        if(mainLoopVar % timestep_denovo == 0)
        {
			addFilamentTouchingNucleatorBead(false);
        }
        
        if(mainLoopVar % timestepGrowFilaments == 0)
        {
            addSegmentsToAllFilaments(numMonomersPerSecondAdd, timeBetweenGrowingFilaments);
        }
        
        if(mainLoopVar % timestepAddBranch == 0)
        {
            addRandomBranch(2);
        }
        
        if(verbose == true)
            cout << "before snapshot" << endl;
        
        
    }
    
    t2=omp_get_wtime();
	double diff ((double)t2-(double)t1);

	cout<<diff<<endl;
    
    outputXYZ.close();
    outputBND.close();
    outputANG.close();
}