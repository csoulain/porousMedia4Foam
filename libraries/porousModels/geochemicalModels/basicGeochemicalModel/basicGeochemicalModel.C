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
      mineralList_(dict.lookup("mineral")),
      Ys_(mineralList_.size() ),
      dMs_(mineralList_.size() ),
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
      porousMedia_(mineralList_.size()),
      absolutePermeabilityModelPtr_
      (
          absolutePermeabilityModel::New(mesh, dict)
      ),
      phiName_(dict.lookupOrDefault<word>("phi","phi")),
      phi_(mesh.lookupObject<surfaceScalarField>(phiName_))
///      rhol_(dict.lookup("rhol"))
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

      dMs_.set
      (
        s,
        new volScalarField
        (
          IOobject
          (
            "dMs."+mineralList_[s],
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ, //READ_IF_PRESENT,  //MUST_READ ??
            IOobject::NO_WRITE
          ),
          mesh_,
          dimensionedScalar(currentMineral,dimless,0.0),
          "zeroGradient"
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
          dict
        )
      );
    }
    updatePorosity();


}

// -------------------------------------------------------------------------//

/*
Foam::volScalarField Foam::basicGeochemicalModel::dMl() const
{

    volScalarField dMl_(0.0*fvc::ddt(Y_[0])/this->rhol());
    forAll(Y_,s)
    {
        dMl_ = dMl_ + fvc::ddt(Y_[s])/this->rhol();
    }

    return dMl_;
}
*/

/*
Foam::volScalarField Foam::basicGeochemicalModel::rhol() const
{

    volScalarField rhol_(0.0*Y_[0]);
    forAll(Y_,s)
    {
        rhol_ = rhol_ + Y_[s];
    }

    return rhol_;
}
*/

void Foam::basicGeochemicalModel::updatePorosity()
{
    eps_ = 0.0*eps_;
    forAll(mineralList_,s)
    {
        eps_+=Ys_[s];
    }
    eps_ = 1.-eps_-inertMineral_;
}



// ************************************************************************* //
