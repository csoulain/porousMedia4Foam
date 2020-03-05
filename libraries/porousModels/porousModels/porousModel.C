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

#include "porousModel.H"
#include "fvcDdt.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porousModel::porousModel
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
        mesh_(mesh),
        porousMediaName_("porousMedia"),
        porousMediaDict_(dict.subDict(porousMediaName_+"Properties")),
        Ys_
        (
            IOobject
            (
                "Ys",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("Ys",dimless,0.0),
            "zeroGradient"
        ),
        Yss_(Ys_),
        eps_
        (
            IOobject
            (
                "eps",
                mesh.time().timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            porousMediaDict_.lookupOrDefault("eps",dimensionedScalar("",dimless,1.))
        ),
        absolutePermeabilityModelPtr_
        (
            absolutePermeabilityModel::New(mesh, porousMediaDict_)
        ),
        surfaceAreaModelPtr_
        (
            surfaceAreaModel::New(mesh, eps_, porousMediaDict_)
        )
{}

Foam::porousModel::porousModel
(
    const fvMesh& mesh,
    const word & name,
    const volScalarField& Ys,
    const dictionary& dict
)
:
        mesh_(mesh),
        porousMediaName_(name),
        porousMediaDict_(dict.subDict(porousMediaName_+"Properties")),
        Ys_(Ys),
        Yss_(Ys),
        eps_
        (
            IOobject
            (
                "eps",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            1.-Ys-SMALL
        ),
        absolutePermeabilityModelPtr_(NULL),
        surfaceAreaModelPtr_
        (
            surfaceAreaModel::New(mesh, Yss_, porousMediaDict_)
        )
{}


// -------------------------------------------------------------------------//


// ************************************************************************* //
