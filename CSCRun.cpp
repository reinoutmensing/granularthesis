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
#include "Mercury3DRestart.h"
#include <iostream>

class CSCRun : public Mercury3DRestart
{
public:

    CSCRun(Mdouble shearVelocity): shearVelocity_(shearVelocity) {};

    void setupInitialConditions() override
    {
        Mdouble timeMax = getTimeMax();
        setName("CSCInit");
        logger(INFO, "Reading file %\n", restartFile.getName(), Flusher::NO_FLUSH);
        readRestartFile();
        setTimeMax(timeMax);
        setTime(0);
        setRestarted(false);
        setName("CSCRun");
        writeXBallsScript();

        //setParticlesWriteVTK(true);// Adds VTU files
        //wallHandler.setWriteVTK(true);// Adds VTU files
        setXBallsAdditionalArguments("-v0 -solidf -cmode 5");
        setFileType(FileType::ONE_FILE);

        logger(INFO, "loaded % fixed particles", particleHandler.getNumberOfObjects());

        // Set velocities for the walls
        for (BaseParticle* p : particleHandler) {
            if (p->isFixed()) {
                if (p->getPosition().X > .0) // Changes the location of the slit
                    p->setVelocity(Vec3D(0.0, 0.5 * shearVelocity_ , 0.0));
                else
                    p->setVelocity(Vec3D(0.0, -0.5 * shearVelocity_, 0.0));
            }
        }
        
        // Save every 50 timesteps
        setSaveCount((0.25 / shearVelocity_) / getTimeStep());
        std::cout << "Time step: " << getTimeStep() << std::endl;
    }

    void printTime() const override
    {
        logger(INFO, "t=% ene=% wallTime=%",
               getTime(), getKineticEnergy() / getElasticEnergy(), getWallTime() - getInitialWallTime());
    }

    Mdouble shearVelocity_;
};

int main(int argc, char *argv[])
{   
    Mdouble shearVelocity =  1.0/40.0; // No shear velocity for a static situation
    CSCRun SC(shearVelocity);
    SC.setSaveCount(100);  
    SC.setMaxWallTime(100 * 3600); // Kill simulation after 30 hours
    SC.setTimeMax(20000.0);        // Maximum simulation time
    SC.solve(argc, argv);
    return 0;
}
