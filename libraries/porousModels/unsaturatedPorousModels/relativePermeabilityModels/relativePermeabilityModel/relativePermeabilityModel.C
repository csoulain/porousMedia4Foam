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

#include "relativePermeabilityModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(relativePermeabilityModel, 0);
defineRunTimeSelectionTable(relativePermeabilityModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::relativePermeabilityModel::relativePermeabilityModel
(
    const word& name,
    const dictionary& dict,
    const volScalarField& Sb,
    const autoPtr<reducedSaturationModel> & reducedSaturationModelPtr
)
    :
    name_(name),
    transportProperties_(dict),
    Sb_(Sb),
    reducedSaturationModelPtr_(reducedSaturationModelPtr),
    Se_(reducedSaturationModelPtr_->Se()),
    dSedS_(reducedSaturationModelPtr_->dSedS()),
    kra_
    (
        IOobject
        (
            name+".kra",
            Sb_.time().timeName(),
            Sb_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Sb.mesh(),
        dimensionSet(0,0,0,0,0)
    ),
    krb_
    (
        IOobject
        (
            name+".krb",
            Sb_.time().timeName(),
            Sb_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Sb.mesh(),
        dimensionSet(0,0,0,0,0)
    ),
    dkradS_
    (
        IOobject
        (
            name+".dkradS",
            Sb_.time().timeName(),
            Sb_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Sb.mesh(),
        dimensionSet(0,0,0,0,0)
    ),
    dkrbdS_
    (
        IOobject
        (
            name+".dkrbdS",
            Sb_.time().timeName(),
            Sb_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Sb.mesh(),
        dimensionSet(0,0,0,0,0)
    ),
    kraf_(fvc::interpolate(kra_,"kra")),
    krbf_(fvc::interpolate(krb_,"krb"))
{}

// ************************************************************************* //
