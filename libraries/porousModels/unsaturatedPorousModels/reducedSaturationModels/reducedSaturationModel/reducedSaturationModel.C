/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "reducedSaturationModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(reducedSaturationModel, 0);
defineRunTimeSelectionTable(reducedSaturationModel, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reducedSaturationModel::reducedSaturationModel
(
    const dictionary& reducedSaturationProperties,
    const volScalarField& Sb
)
    :
    reducedSaturationProperties_(reducedSaturationProperties),
    Sb_(Sb),
    Se_
    (
        IOobject
        (
            "Se",
            Sb.time().timeName(),
            Sb.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Sb,
        calculatedFvPatchScalarField::typeName
    ),
    dSedS_
    (
        IOobject
        (
            "dSedS",
            Sb.time().timeName(),
            Sb.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Sb.mesh(),
        dimensionedScalar("dSedS",dimless,1.0),
        calculatedFvPatchScalarField::typeName
    )
{}

// ************************************************************************* //
