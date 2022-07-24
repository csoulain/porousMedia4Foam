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

#include "minMax.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace reducedSaturationModels
{
defineTypeNameAndDebug(minMax, 0);

addToRunTimeSelectionTable
(
    reducedSaturationModel,
    minMax,
    dictionary
);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reducedSaturationModels::minMax::minMax
(
    const dictionary& dict,
    const volScalarField& Sb
)
    :
    reducedSaturationModel(dict,Sb),
    minMaxCoeffs_(dict.subDict(typeName + "Coeffs")),
    Smin_
    (
        IOobject
        (
            Sb_.name()+"min",
            Sb_.time().timeName(),
            Sb_.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        Sb.mesh(),
        minMaxCoeffs_.lookupOrDefault(Sb_.name()+"min",dimensionedScalar(Sb_.name()+"min",dimless,0))
    ),
    Smax_
    (
        IOobject
        (
            Sb_.name()+"max",
            Sb_.time().timeName(),
            Sb_.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        Sb.mesh(),
        minMaxCoeffs_.lookupOrDefault(Sb_.name()+"max",dimensionedScalar(Sb_.name()+"max",dimless,0))
    )
{
    Se_ = (Sb_-Smin_)/(Smax_-Smin_);

    Info << "Parameters for reduced saturation model" << nl << "{" << endl;
    Info << "    Smax ";
    if (Smax_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(Smax_).value() << endl;}
    Info <<  "    Smin ";
    if (Smin_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(Smin_).value() << endl;}
    Info << "} \n" << endl;
}

// ************************************************************************* //
