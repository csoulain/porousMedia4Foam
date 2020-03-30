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

#include "mingliangLaw.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceAreaModels
{
    defineTypeNameAndDebug(mingliangLaw, 0);

    addToRunTimeSelectionTable
    (
        surfaceAreaModel,
        mingliangLaw,
        dictionary
    );
}
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceAreaModels::mingliangLaw::mingliangLaw
(
    const fvMesh& mesh,
    const volScalarField& Ys,
    const dictionary& dict
)
:
    surfaceAreaModel(mesh, Ys, dict),
    mingliangLawDict_(dict.subDict(typeName+"Coeffs")),
    A0_(mingliangLawDict_.lookupOrDefault
    (
        "A0",
        dimensionedScalar("A0",dimensionSet(0,-1,0,0,0,0,0),1.0/*0.0*/))
    ),
    n_(readScalar(mingliangLawDict_.lookup("n"))),
    initVolFrac_(readScalar(mingliangLawDict_.lookup("initVF"))),
    Ae_
    (
        IOobject
        (
          "Ae",
          mesh_.time().timeName(),
          mesh_,
          IOobject::READ_IF_PRESENT,
          IOobject::AUTO_WRITE
        ),
        mesh_,
        A0_,
        "zeroGradient"
    ),
    Ys_(Ys)
{}

// * * * * * * * * * * * * * * member functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::surfaceAreaModels::mingliangLaw::surfaceArea() const
{
      return Ae_;
}

void Foam::surfaceAreaModels::mingliangLaw::updateSurfaceArea()
{
    if (scalar(initVolFrac_) == 0.)
     Ae_ = A0_;
    else
    {
     //Info << "Mineral VF =" << Ys_ << endl;
     Ae_ = A0_ * Foam::pow(Ys_/scalar(initVolFrac_),scalar(n_));
     //Info << "Ae = " << Ae_ << nl;
    }
}

// -------------------------------------------------------------------------//
