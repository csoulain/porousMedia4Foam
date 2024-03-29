/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    darcyFoam

Description
    Stationary solver for incompressible single-phase flow in porous medium

Developers
    - Pierre Horgue
    - 10/02/2020 - CS : modified with the porousModel class
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "incompressiblePhase.H"
#include "geochemicalModel.H"
//#include "porousModel.H"
#include "sourceEventFile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

//    argList::addOption("phase","a","specify the phase name");
//    Foam::argList args(argc,argv);
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readGravitationalAcceleration.H"
    #include "createFields.H"
    #include "readEvent.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {

        Info << "Time = " << runTime.timeName() << nl << endl;

        #include "CourantNo.H"

        porousMedia.update();


        Mf = Kf / mu;

        fvScalarMatrix pEqn
        (
            fvm::laplacian(-Mf,p) + fvc::div(phiG) - sourceTerm
        );

        pEqn.solve();

        phi = pEqn.flux() + phiG;

        U = fvc::reconstruct(phi);
        U.correctBoundaryConditions();

        runTime.write();

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
