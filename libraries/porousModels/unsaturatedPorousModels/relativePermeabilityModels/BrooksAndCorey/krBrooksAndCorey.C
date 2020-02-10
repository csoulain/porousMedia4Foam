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

#include "krBrooksAndCorey.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace relativePermeabilityModels
{
defineTypeNameAndDebug(krBrooksAndCorey, 0);

addToRunTimeSelectionTable
(
    relativePermeabilityModel,
    krBrooksAndCorey,
    dictionary
);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::relativePermeabilityModels::krBrooksAndCorey::krBrooksAndCorey
(
    const word& name,
    const dictionary& dict,
    const volScalarField& Sb,
    const autoPtr<reducedSaturationModel> & reducedSaturationModelPtr
)
    :
    relativePermeabilityModel(name, dict,Sb,reducedSaturationModelPtr),
//    Se_(Sb),
    krBrooksAndCoreyCoeffs_(dict.subDict(typeName + "Coeffs")),
    n_
    (
        IOobject
        (
            "n",
            Se_.time().timeName(),
            Se_.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        Sb.mesh(),
        dimensionedScalar("n",dimless,krBrooksAndCoreyCoeffs_.lookupOrDefault<scalar>("n",0))
    ),
    kramax_
    (
        IOobject
        (
            "kr"+Se_.name()+"max",
            Se_.time().timeName(),
            Se_.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        Sb.mesh(),
        dimensionedScalar("kr"+Se_.name()+"max",dimless,krBrooksAndCoreyCoeffs_.lookupOrDefault<scalar>("kr"+Se_.name()+"max",1.0))
    ),
    krbmax_
    (
        IOobject
        (
            "kr"+Se_.name()+"max",
            Se_.time().timeName(),
            Se_.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        Sb.mesh(),
        dimensionedScalar("kr"+Se_.name()+"max",dimless,krBrooksAndCoreyCoeffs_.lookupOrDefault<scalar>("kr"+Se_.name()+"max",1.0))
    )
{

    if (gMin(n_) <= 0)
    {
        FatalErrorIn
            (
                "in krBrooksAndCorey.C"
            )
            << "Relative permeability coefficient n equal or less than 0"
                << exit(FatalError);
    }

    Info << "Brooks and Corey parameters for relative permeability model" << nl << "{" << endl;
    Info << "    n ";
    if (n_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(n_).value() << endl;}
    Info << "    kramax ";
    if (kramax_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(kramax_).value() << endl;}
    Info << "    krbmax ";
    if (krbmax_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(krbmax_).value() << endl;}
    Info << "} \n" << endl;
}

// ************************************************************************* //
