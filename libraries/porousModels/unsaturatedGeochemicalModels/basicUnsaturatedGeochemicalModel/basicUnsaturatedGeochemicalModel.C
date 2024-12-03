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

#include "basicUnsaturatedGeochemicalModel.H"
#include "fvcDdt.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(basicUnsaturatedGeochemicalModel, 0);
    defineRunTimeSelectionTable(basicUnsaturatedGeochemicalModel, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::basicUnsaturatedGeochemicalModel::basicUnsaturatedGeochemicalModel
(
    const fvMesh& mesh,
    const dictionary& dict,
    const volScalarField &Sb,
    const incompressiblePhase &phasea,
    const incompressiblePhase &phaseb
)
:
      fluidProperties(mesh,dict),
      unsaturatedPorousAssemblage(mesh, dict, Sb, phasea, phaseb),
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
      dGasvdRho_
      (
          IOobject
          (
              "dGasvdRho",
              mesh.time().timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::NO_WRITE
          ),
          mesh,
          dimensionedScalar("dGasvdRho",dimless/dimTime,0.0),
          "zeroGradient"
      ),
      phiaName_(dict.lookupOrDefault<word>("phia","phia")),
      phibName_(dict.lookupOrDefault<word>("phib","phib")),
      phia_(mesh.lookupObject<surfaceScalarField>(phiaName_)),
      phib_(mesh.lookupObject<surfaceScalarField>(phibName_)),
      phiName_(geochemicalModelDict_.lookupOrDefault<word>("phi","phi")),
//      phi_(mesh.lookupObject<surfaceScalarField>(phiName_)),
      phi_(phia_+phib_),
      SbName_(dict.lookupOrDefault<word>("Sb","Sb")),
      Sb_(mesh.lookupObject<volScalarField>(SbName_))

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

void Foam::basicUnsaturatedGeochemicalModel::updatedMinvdRho()
{
    dMinvdRho_ = 0.0*dMinvdRho_;
 //   forAll(mineralList_,s)
    {
        //dMinvdRho_+= -rhos_[s]*fvc::ddt(Ys_[s])*(1./rhol_-1./rhos_[s]);
//        dMinvdRho_+= -rhos_[s]*fvc::ddt(Ys_[s])*(1./this->rho()-1./rhos_[s]);
    }
    dMinvdRho_.correctBoundaryConditions(); //necessary??
}

void Foam::basicUnsaturatedGeochemicalModel::updatedGasvdRho()
{
    dGasvdRho_ = 0.0*dGasvdRho_;
 //   forAll(mineralList_,s)
    {
        //dMinvdRho_+= -rhos_[s]*fvc::ddt(Ys_[s])*(1./rhol_-1./rhos_[s]);
//        dMinvdRho_+= -rhos_[s]*fvc::ddt(Ys_[s])*(1./this->rho()-1./rhos_[s]);
    }
    dGasvdRho_.correctBoundaryConditions(); //necessary??
}


// ************************************************************************* //
