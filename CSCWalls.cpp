//Copyright (c) 2013-2023, The MercuryDPM Developers Team. All rights reserved.
//For the list of developers, see <http://www.MercuryDPM.org/Team>.
//
//Redistribution and use in source and binary forms, with or without
//modification, are permitted provided that the following conditions are met:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name MercuryDPM nor the
//    names of its contributors may be used to endorse or promote products
//    derived from this software without specific prior written permission.
//
//THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//DISCLAIMED. IN NO EVENT SHALL THE MERCURYDPM DEVELOPERS TEAM BE LIABLE FOR ANY
//DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
//LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
//ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

//based on /storage2/usr/people/sluding/MDCC/C3DshearXL30/MU0_LONG2
#include "Mercury3D.h"
#include "Species/LinearViscoelasticSlidingFrictionSpecies.h"
#include "Boundaries/PeriodicBoundary.h"
#include "Walls/AxisymmetricIntersectionOfWalls.h"
#include <iostream>

class CSCWalls : public Mercury3D {
public:
    CSCWalls (Mdouble width, Mdouble length, Mdouble height, Mdouble wallThickness)
    {
        eneTolerance = 1.0;
        sizeDistribution = cbrt(2.0);
    
        setXMin(-width / 2 - wallThickness);
        setXMax(width / 2 + wallThickness);
        setYMin(0);
        setYMax(length);
        setZMin(-wallThickness);
        setZMax(height);
    
        logger(INFO, "Creating flat-wall box of"
                     " width [%, %]"
                     " length [%, %]"
                     " height [%, %]\n", Flusher::NO_FLUSH, getXMin(), getXMax(), getYMin(), getYMax(), getZMin(),
               getZMax());
    
        //set default name
        setName("CSCWalls");
        setXBallsAdditionalArguments("-v0 -solidf -3dturn 0");
    
        //set gravity (extra strong since we have few particle layers)
        setGravity(Vec3D(0.0, 0.0, -1));
    
        //create new species
        auto particleSpecies = speciesHandler.copyAndAddObject(LinearViscoelasticSlidingFrictionSpecies());
        particleSpecies->setDensity(6.0 / constants::pi);
    
        //set inter-particle contact properties
        particleSpecies->setStiffness(10e5 / 100.0);
        particleSpecies->setDissipation(25.0 / 10.0);
        particleSpecies->setSlidingStiffness(2.0 / 7.0 * particleSpecies->getStiffness());
        particleSpecies->setSlidingDissipation(particleSpecies->getDissipation());
        particleSpecies->setSlidingFrictionCoefficient(0.5);

        //create new (particle) species
        wallSpecies = speciesHandler.copyAndAddObject(particleSpecies);
        particleSpecies->setSlidingFrictionCoefficient(0.2);
        auto particleWallSpecies = speciesHandler.getMixedObject(particleSpecies, wallSpecies);
        particleWallSpecies->setSlidingFrictionCoefficient(constants::inf);
        
        setTimeStep(0.02 * particleSpecies->getCollisionTime(particleSpecies->getMassFromRadius(0.5)));
        logger(INFO, "recommended time step: %\n", getTimeStep(), Flusher::NO_FLUSH);
        setTimeStep(5.0 * getTimeStep());
        logger(INFO, "time step use for creating initial conditions: %", getTimeStep());

        //particleSpecies->getCollisionTime(particleSpecies->getMassFromRadius(0.5))
        //particleSpecies->getMassFromRadius(0.5)
        //radius = 0.5
        //slidingfriction = inf
        //
        // Calculate mass from radius 0.5
        Mdouble mass = particleSpecies->getMassFromRadius(0.5);
         // Retrieve and print predefined constants
        // Print values using the correct function calls and passing mass as needed
        std::cout << "Restitution: " << particleSpecies->getRestitutionCoefficient(mass) << std::endl;
        std::cout << "Sliding Friction: " << particleSpecies->getSlidingFrictionCoefficient() << std::endl;
        std::cout << "Body Force: " << getGravity().Z << std::endl;
        std::cout << "Mass for radius 0.5: " << mass << std::endl;
        std::cout << "Collision Time for mass 0.5: "
                << particleSpecies->getCollisionTime(mass) << std::endl;
    
        //set time step
        setTimeMax(1e20);
        setSaveCount(2000);
    
        //set boundaries
        PeriodicBoundary b;
        b.set(Vec3D(0, 1, 0), getYMin(), getYMax());
        boundaryHandler.copyAndAddObject(b);

        // And in order to make periodic boundary of x this should be included.
        // b.set(Vec3D(1,0,0), getXMin(), getXMax());
        // boundaryHandler.copyAndAddObject(b);

        //set walls
        InfiniteWall w;
        w.setSpecies(wallSpecies);

        // in order to make periodic boundary of x this can be commented out
        w.set(Vec3D(1,0,0), getMax());
        wallHandler.copyAndAddObject(w);
        w.set(Vec3D(-1,0,0), getMin());
        wallHandler.copyAndAddObject(w);


        //try to make the sides into an infinitewall as well. 
        w.set(Vec3D(0,0,-1), getMin());
        wallHandler.copyAndAddObject(w);

        
        

        //remove wall thickness
        setXMin(-width/2);
        setXMax( width/2);
        setZMin(0);
    }

    // Getter function to access particleSpecies
    LinearViscoelasticSlidingFrictionSpecies* getParticleSpecies() const {
        return wallSpecies; // Or particleSpecies if that's the correct variable
    }


    bool continueSolve() const override
    {
        static unsigned int counter = 0;
        if (++counter>100)
        {
            counter=0;
            if (getKineticEnergy()<eneTolerance*getElasticEnergy())
                return false;
        }
        return true;
    }

    void printTime() const override
    {
        logger(INFO, "t=% Ene=%", getTime(), getKineticEnergy() / getElasticEnergy());
    }

    //add flow particles
    void setupInitialConditions() override
    {
        //number of particles to be inserted
        unsigned int n = (getXMax() - getXMin()) * (getYMax() - getYMin())
                         * 1.8 * (getZMax() - getZMin());
        logger(INFO, "Inserting % particles\n", n, Flusher::NO_FLUSH);
        //try to find new insertable particles
        unsigned int i = 0;
        SphericalParticle p;
        p.setSpecies(wallSpecies);
        Mdouble s = sizeDistribution;
        Mdouble rMin = cbrt(0.5 / (s * s + 1.0) / (s + 1.0));
        p.setRadius(s * rMin);
        logger(INFO, "Particle sizes from % to %  (sizeDistribution %)\n",
               Flusher::NO_FLUSH, rMin, s * rMin, sizeDistribution);
        Vec3D position;
        Mdouble a = 0.0;
        while (i < n)
        {
            position.X = random.getRandomNumber(getXMin() + p.getRadius(),
                                                getXMax() - p.getRadius());
            position.Y = random.getRandomNumber(getYMin() + p.getRadius(),
                                                getYMax() - p.getRadius());
            position.Z = random.getRandomNumber(getZMin() + p.getRadius(),
                                                std::max(getZMin() + p.getRadius(), a * getZMax() - p.getRadius()));
            p.setPosition(position);
            if (checkParticleForInteraction(p))
            {
                //logger(INFO, "%", i);
                particleHandler.copyAndAddObject(p);
                p.setRadius(random.getRandomNumber(rMin, s * rMin));
                i++;
                if (particleHandler.getNumberOfObjects() % 100 == 0)
                {
                    logger(INFO, " %", particleHandler.getNumberOfObjects() / 100);
                }
            
            }
            else
            {
                a += 0.00001;
            }
        }
        logger(INFO, "\nInserted % particles within z<%", n, a * getZMax());
    }

    void saveWalls()
    {
        logger(INFO, "Creating flat-wall box of"
                     " width [%, %]"
                     " length [%, %]"
                     " height [%, %]\n",
               Flusher::NO_FLUSH, getXMin(), getXMax(), getYMin(), getYMax(), getZMin(), getZMax());
    
        //only keep wall particles
        for (BaseParticle* p: particleHandler)
        {
            Mdouble x = p->getPosition().X;
            Mdouble z = p->getPosition().Z;
            if (x > getXMax() || x < getXMin() || z < getZMin())
            {
                if (z < getZMax())
                {
                    p->fixParticle();
                }
            }
        }
        for (int i = particleHandler.getNumberOfObjects() - 1; i >= 0; i--)
        {
            if (!particleHandler.getObject(i)->isFixed())
                particleHandler.removeObject(i);
        }
        interactionHandler.clear();
    
        logger(INFO, "kept % fixed particles", particleHandler.getNumberOfObjects());
    
        //save data and restart file of last time step
        dataFile.open();
        outputXBallsData(dataFile.getFstream());
        dataFile.close();
    
        restartFile.open();
        writeRestartFile();
        restartFile.close();
    }

    LinearViscoelasticSlidingFrictionSpecies* wallSpecies;
    Mdouble eneTolerance;
    Mdouble sizeDistribution;
};

int main(int argc, char *argv[]) {
    Mdouble width = 30; //in particle diameter
    Mdouble length = 10;//in particle diameter
    Mdouble height = 150;//in particle diameter
    Mdouble wallThickness = 1;//in particle diameter
    CSCWalls SC(width, length, height, wallThickness);

   
    SC.solve(argc, argv);
    SC.saveWalls();
    return 0;
}
