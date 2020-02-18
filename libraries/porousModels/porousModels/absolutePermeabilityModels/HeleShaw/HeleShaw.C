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

#include "HeleShaw.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace absolutePermeabilityModels
{
    defineTypeNameAndDebug(HeleShaw, 0);

    addToRunTimeSelectionTable
    (
        absolutePermeabilityModel,
        HeleShaw,
        dictionary
    );
}
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::absolutePermeabilityModels::HeleShaw::HeleShaw
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    absolutePermeabilityModel(mesh, dict),
    HeleShawDict_(dict.subDict(typeName+"Coeffs")),
    h_(HeleShawDict_.lookup("h")),
    K_
    (
        IOobject
        (
          "K",
          mesh_.time().timeName(),
          mesh_,
          IOobject::NO_READ,
          IOobject::NO_WRITE
        ),
        mesh_,
        pow(h_,2)/12,
        "zeroGradient"
    ),
    invK_
    (
        IOobject
        (
          "invK",
          mesh_.time().timeName(),
          mesh_,
          IOobject::NO_READ,
          IOobject::NO_WRITE
        ),
        1/K_,
        "zeroGradient"
    ),
    Kf_("Kf", fvc::interpolate(K_,"K"))
{}

// * * * * * * * * * * * * * * member functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::absolutePermeabilityModels::HeleShaw::absolutePermeability() const
{
      return K_;
}

Foam::tmp<Foam::volScalarField>
Foam::absolutePermeabilityModels::HeleShaw::inversePermeability() const
{
      return invK_;
}

Foam::tmp<Foam::surfaceScalarField>
Foam::absolutePermeabilityModels::HeleShaw::Kf() const
{
      return Kf_;
}

void Foam::absolutePermeabilityModels::HeleShaw::updatePermeability()
{
    //do nothing

}

// -------------------------------------------------------------------------//
