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
      mesh_(mesh),
      geochemicalModelDict_(dict.subDict("geochemicalProperties")),
      fluidPropertiesDict_(dict.subDict("fluidProperties")),
      mineralList_(geochemicalModelDict_.lookup("mineral")),
      Ys_(mineralList_.size() ),
      inertMineral_
      (
          IOobject
          (
              "inertMineral",
              mesh.time().timeName(),
              mesh,
              IOobject::READ_IF_PRESENT,
              IOobject::AUTO_WRITE
          ),
          mesh,
          dimensionedScalar("inertMineral",dimless,0.0),
          "zeroGradient"
      ),
      eps_
      (
          IOobject
          (
              "eps",
              mesh.time().timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::AUTO_WRITE
          ),
          mesh,
          dimensionedScalar("eps",dimless,1.0),
          "zeroGradient"
      ),
      eps0_
      (
          IOobject
          (
              "eps0",
              mesh.time().timeName(),
              mesh,
              IOobject::READ_IF_PRESENT,
              IOobject::NO_WRITE
          ),
          mesh,
          dimensionedScalar("eps0",dimless,1.0),
          "zeroGradient"
      ),
      rhos_(mineralList_.size() ),
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
      porousMedia_(mineralList_.size()),
      densityModelPtr_
      (
          densityModel::New(mesh, fluidPropertiesDict_)
      ),
  //    rho_(this->rho()),
      /*
      rho_
      (
          IOobject
          (
              "rho",
              mesh.time().timeName(),
              mesh,
              IOobject::READ_IF_PRESENT,
              IOobject::NO_WRITE
          ),
          mesh,
          dimensionedScalar("rho",dimensionSet(1,-3,0,0,0,0,0),1.e3),
          "zeroGradient"
      ),
      */
      absolutePermeabilityModelPtr_
      (
          absolutePermeabilityModel::New(mesh, geochemicalModelDict_)
      ),

      dispersionTensorModelPtr_
      (
          dispersionTensorModel::New(mesh, geochemicalModelDict_)
      ),
      phiName_(geochemicalModelDict_.lookupOrDefault<word>("phi","phi")),
      phi_(mesh.lookupObject<surfaceScalarField>(phiName_))
{

    forAll(mineralList_,s)
    {
      word currentMineral = mineralList_[s];
      Info << " Doing stuff for mineral: " << currentMineral << endl;

      Ys_.set
      (
        s,
        new volScalarField
        (
          IOobject
          (
            "Ys."+mineralList_[s],
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ, //READ_IF_PRESENT,  //MUST_READ ??
            IOobject::AUTO_WRITE
          ),
          mesh_ //,
          //		dimensionedScalar(currentMineral,dimless,0.0),
          //		"zeroGradient"
        )
      );
      Ys_[s].write();

      rhos_.set
      (
          s,
          new dimensionedScalar
          (
              geochemicalModelDict_.subDict(currentMineral+"Properties").lookup("rhos")
          )
      );

      porousMedia_.set
      (
        s,
        new porousModel
        (
          mesh,
          mineralList_[s],
          Ys_[s],
          geochemicalModelDict_
        )
      );
    }
    updatePorosity();


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
}


void Foam::basicGeochemicalModel::updatedMinvdRho()
{
    dMinvdRho_ = 0.0*dMinvdRho_;
    forAll(mineralList_,s)
    {
        //dMinvdRho_+= -rhos_[s]*fvc::ddt(Ys_[s])*(1./rhol_-1./rhos_[s]);
        dMinvdRho_+= -rhos_[s]*fvc::ddt(Ys_[s])*(1./this->rho()-1./rhos_[s]);
    }
}



// ************************************************************************* //
