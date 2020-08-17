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

#include "volumeOfSolid.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceAreaModels
{
    defineTypeNameAndDebug(volumeOfSolid, 0);

    addToRunTimeSelectionTable
    (
        surfaceAreaModel,
        volumeOfSolid,
        dictionary
    );
}
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceAreaModels::volumeOfSolid::volumeOfSolid
(
    const fvMesh& mesh,
    const volScalarField& Ys,
    const dictionary& dict
)
:
    surfaceAreaModel(mesh, Ys, dict),
    Ys_(Ys),
    Ae_
    (
        IOobject
        (
          "Ae",
          mesh_.time().timeName(),
          mesh_,
          IOobject::NO_READ,
          IOobject::NO_WRITE
        ),
        mag(fvc::grad(Ys_))
    )
{}

// * * * * * * * * * * * * * * member functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::surfaceAreaModels::volumeOfSolid::surfaceArea() const
{
      return Ae_;
}

void Foam::surfaceAreaModels::volumeOfSolid::updateSurfaceArea()
{
      Ae_= mag(fvc::grad(Ys_));
//      Ae_=Ae_*2.*Ys_;  //(Diffuse interface function)
      Ae_=Ae_*2.*Foam::pow(Ys_,1);  //(Diffuse interface function)

//      Ae_ = 2.0*fvc::average( fvc::interpolate(Ae_, "Ae") );

}

// -------------------------------------------------------------------------//
