/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "basicGeochemicalModel.H"
#include "fvcDdt.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(basicGeochemicalModel, 0);
    defineRunTimeSelectionTable(basicGeochemicalModel, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::basicGeochemicalModel::basicGeochemicalModel
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
      fluidProperties(mesh,dict),
      porousAssemblage(mesh,dict.subDict("geochemicalProperties")),
      mesh_(mesh),
      geochemicalModelDict_(dict.subDict("geochemicalProperties")),
      fluidPropertiesDict_(dict.subDict("fluidProperties")),
      dMinvdRho_
      (
          IOobject
          (
              "dMinvdRho",
              mesh.time().timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::NO_WRITE
          ),
          mesh,
          dimensionedScalar("dMinvdRho",dimless/dimTime,0.0),
          "zeroGradient"
      ),
      phiName_(geochemicalModelDict_.lookupOrDefault<word>("phi","phi")),
      phi_(mesh.lookupObject<surfaceScalarField>(phiName_))
{


// -----------------------------------------------------------------------------
    const word densitymodelType
    (
      fluidPropertiesDict_.lookup("densityModel")
    );

    const word geochemicalmodelType
    (
      dict.lookup("geochemicalModel")
    );

    if
    (
      (densitymodelType == "fromPhreeqc") && (geochemicalmodelType != "phreeqcRM")
    )
    {
      FatalErrorInFunction
          << "fromPhreeqc densityModel type must be used "
          << "with phreeqcRM geochemical package" <<nl
          << exit(FatalError);
    }

// -----------------------------------------------------------------------------


}

// -------------------------------------------------------------------------//


void Foam::basicGeochemicalModel::updatePorosity()
{
    eps_ = 0.0*eps_;
    forAll(mineralList_,s)
    {
        eps_+=Ys_[s];
    }
    eps_ = 1.-eps_-inertMineral_;
    eps_.correctBoundaryConditions(); //necessary??
}


void Foam::basicGeochemicalModel::updatedMinvdRho()
{
    dMinvdRho_ = 0.0*dMinvdRho_;
    forAll(mineralList_,s)
    {
        //dMinvdRho_+= -rhos_[s]*fvc::ddt(Ys_[s])*(1./rhol_-1./rhos_[s]);
        dMinvdRho_+= -rhos_[s]*fvc::ddt(Ys_[s])*(1./this->rho()-1./rhos_[s]);
    }
    dMinvdRho_.correctBoundaryConditions(); //necessary??
}



// ************************************************************************* //
