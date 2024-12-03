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
    impesFoam

Description
    Transient solver for incompressible two-phase flow (Darcy's law) in porous media
    using the IMPES method (IMplicit Pressure Explicit Saturation).
    Permeability is isotropic (K == volScalarField)

Developers
    P. Horgue, C. Soulaine, J. Franc, R. Guibert and G. Debenest
    "An open-source toolbox for multiphase flow in porous media"

    - 2012 : CS - first impesFoam version
    - 2013 : PH - impesFoam part of porousMultiphaseFoam Toolbox
    - 2015 : JF - Coats and Todd time conditions
    - 02/02/2020 - CS: unsaturatedPorousModel class, part of porousMedia4Foam

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "harmonic.H"
#include "incompressiblePhase.H"
#include "unsaturatedGeochemicalModel.H"


#include "uniformDimensionedFields.H"

#include "fvcSnGrad.H"
#include "fvcFlux.H"
#include "fvcSurfaceIntegrate.H"
#include "fvcReconstruct.H"

#include "fvmDdt.H"
//#include "fvmDiv.H"
#include "fvmLaplacian.H"

#include "fixedValueFvPatchField.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createTimeControls.H"
    #include "readGravitationalAcceleration.H"
    #include "createFields.H"
    #include "readTimeControls.H"
//    #include "readEvent.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
  //      if (outputEventIsPresent) outputEvent.updateIndex(runTime.timeOutputValue());
  //      if (sourceEventIsPresent) sourceEvent.updateIndex(runTime.timeOutputValue());
  //      forAll(patchEventList,patchEventi) patchEventList[patchEventi]->updateIndex(runTime.timeOutputValue());
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

  //      #include "computeSourceTerm.H"

        //- Solve saturation equation (explicit)
        #include "SEqn.H"

        porousMedia.update();

        //- Solve pressure equation (implicit)
        #include "pEqn.H"

  //      #include "eventWrite.H"

        if(runTime.writeTime())
        {
          runTime.write();
        }


        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
