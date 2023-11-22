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

#include "constantViscosity.H"
#include "addToRunTimeSelectionTable.H"


//#include "fvCFD.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(constantViscosity, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        constantViscosity,
        dictionary
    );
}
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::constantViscosity::constantViscosity
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    viscosityModel(mesh, dict),
    constantViscosityDict_(dict.subDict(typeName+"Coeffs")),
    nu0_(constantViscosityDict_.lookup("nu0")),
    nu_
    (
        IOobject
        (
          "nu",
          mesh_.time().timeName(),
          mesh_,
          IOobject::READ_IF_PRESENT,
          IOobject::NO_WRITE
        ),
        mesh_,
        nu0_,
        "zeroGradient"
    )
{}

// * * * * * * * * * * * * * * member functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::constantViscosity::nu() const
{
      return nu_;
}

void Foam::viscosityModels::constantViscosity::updateViscosity()
{
    //do nothing

}

// -------------------------------------------------------------------------//
