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

#include "diffusionOnlyTensor.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dispersionTensorModels
{
    defineTypeNameAndDebug(diffusionOnlyTensor, 0);

    addToRunTimeSelectionTable
    (
        dispersionTensorModel,
        diffusionOnlyTensor,
        dictionary
    );
}
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dispersionTensorModels::diffusionOnlyTensor::diffusionOnlyTensor
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    dispersionTensorModel(mesh, dict),
    diffusionOnlyTensorDict_(dict.subDict(typeName+"Coeffs")),
    Di_( diffusionOnlyTensorDict_.lookup("Di") ),
    Deff_
     (
      IOobject
      (
       "Deff",
       mesh_.time().timeName(),
       mesh_,
       IOobject::READ_IF_PRESENT,
       IOobject::NO_WRITE
      ),
      mesh_,
      dimensionedTensor("Deff",dimensionSet(0, 2, -1, 0, 0),tensor::zero),
      "zeroGradient"
    )
{
}

// * * * * * * * * * * * * * * member functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volTensorField>
Foam::dispersionTensorModels::diffusionOnlyTensor::effectiveDispersionTensor() const
{
      return Deff_;
}

void Foam::dispersionTensorModels::diffusionOnlyTensor::updateDispersionTensor()
{
       Deff_= tensor(1,0,0,0,1,0,0,0,0)*Di_;
}
// -------------------------------------------------------------------------//
