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

#include "heterogeneousScalarConstant.H"
#include "addToRunTimeSelectionTable.H"


//#include "fvCFD.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace absolutePermeabilityModels
{
    defineTypeNameAndDebug(heterogeneousScalarConstant, 0);

    addToRunTimeSelectionTable
    (
        absolutePermeabilityModel,
        heterogeneousScalarConstant,
        dictionary
    );
}
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::absolutePermeabilityModels::heterogeneousScalarConstant::heterogeneousScalarConstant
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    absolutePermeabilityModel(mesh, dict),
    heterogeneousScalarConstantDict_(dict.subDict(typeName+"Coeffs")),
    K0_(heterogeneousScalarConstantDict_.lookupOrDefault
    (
        "K0",
        dimensionedScalar("K0",dimensionSet(0,2,0,0,0,0,0),0))
    ),
    K_
    (
        IOobject
        (
          "K",
          mesh_.time().constant(),
          mesh_,
          IOobject::READ_IF_PRESENT,
          IOobject::NO_WRITE
        ),
        mesh_,
        K0_,
        "zeroGradient"
    )
{

}

// * * * * * * * * * * * * * * member functions  * * * * * * * * * * * * * * //
/*
Foam::tmp<Foam::volScalarField>
Foam::absolutePermeabilityModels::heterogeneousScalarConstant::inversePermeability() const
{
      return invK_;
}
*/
Foam::tmp<Foam::volScalarField>
Foam::absolutePermeabilityModels::heterogeneousScalarConstant::absolutePermeability() const
{
      return K_;
}

//void Foam::absolutePermeabilityModels::heterogeneousScalarConstant::updatePermeability()
//{}

// -------------------------------------------------------------------------//
