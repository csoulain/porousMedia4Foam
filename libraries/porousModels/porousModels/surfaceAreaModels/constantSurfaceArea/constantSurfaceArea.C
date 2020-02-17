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

#include "constantSurfaceArea.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceAreaModels
{
    defineTypeNameAndDebug(constantSurfaceArea, 0);

    addToRunTimeSelectionTable
    (
        surfaceAreaModel,
        constantSurfaceArea,
        dictionary
    );
}
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceAreaModels::constantSurfaceArea::constantSurfaceArea
(
    const fvMesh& mesh,
    const volScalarField& eps,
    const dictionary& dict
)
:
    surfaceAreaModel(mesh, eps, dict),
    constantSurfaceAreaDict_(dict.subDict(typeName+"Coeffs")),
    A0_(constantSurfaceAreaDict_.lookupOrDefault
    (
        "A0",
        dimensionedScalar("A0",dimensionSet(0,-1,0,0,0,0,0),0.0))
    ),
    Ae_
    (
        IOobject
        (
          "Ae",
          mesh_.time().timeName(),
          mesh_,
          IOobject::READ_IF_PRESENT,
          IOobject::NO_WRITE
        ),
        mesh_,
        A0_,
        "zeroGradient"
    )
{}

// * * * * * * * * * * * * * * member functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::surfaceAreaModels::constantSurfaceArea::surfaceArea() const
{
      return Ae_;
}

void Foam::surfaceAreaModels::constantSurfaceArea::updateSurfaceArea()
{
    //do nothing

}

// -------------------------------------------------------------------------//
