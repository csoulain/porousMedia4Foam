/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    PDRFoam

Description


\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"


#include "pimpleControl.H"
#include "pressureReference.H"
#include "findRefCell.H"
#include "constrainPressure.H"
#include "constrainHbyA.H"
#include "adjustPhi.H"
#include "uniformDimensionedFields.H"
#include "fvModels.H"
#include "fvConstraints.H"

#include "fvcDdt.H"
#include "fvcGrad.H"
#include "fvcFlux.H"
#include "fvcReconstruct.H"
#include "fvcMeshPhi.H"

#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"


#include "fvcSmooth.H"
#include "geochemicalModel.H"
#include "HeleShaw.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "readGravitationalAcceleration.H"
    #include "createFields.H"
   // #include "createFieldRefs.H"
    #include "initContinuityErrs.H"
    #include "createTimeControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    //turbulence->validate();
    scalar StCoNum = 0.0;

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (pimple.run(runTime))
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        runTime++;
        Info<< "\n\nTime = " << runTime.name() << endl;

//        #include "rhoEqn.H"

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            fvModels.correct();

            if (pimple.predictTransport())
            {
     //           turbulence->predict();
     //           thermophysicalTransport.predict();
            }

            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.correctTransport())
            {
  //              turbulence->correct();
  //              thermophysicalTransport.correct();
                  geochemistry.update();
            }
        }

        runTime.write();

        Info<< "\nExecutionTime = "
             << runTime.elapsedCpuTime()
             << " s\n" << endl;
    }

    Info<< "\n end\n";

    return 0;
}


// ************************************************************************* //
