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

#include "basicGeochemicalModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(basicGeochemicalModel, 0);
    defineRunTimeSelectionTable(basicGeochemicalModel, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::basicGeochemicalModel::basicGeochemicalModel
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
//      basicGeochemicalModel(mesh, dict.subDict("geochemicalProperties")),
      mesh_(mesh),
      mineralList_(dict.lookup("mineral"))
///      rhol_(dict.lookup("rhol"))
{
}

// -------------------------------------------------------------------------//

/*
Foam::volScalarField Foam::basicGeochemicalModel::dMl() const
{

    volScalarField dMl_(0.0*fvc::ddt(Y_[0])/this->rhol());
    forAll(Y_,s)
    {
        dMl_ = dMl_ + fvc::ddt(Y_[s])/this->rhol();
    }

    return dMl_;
}
*/

/*
Foam::volScalarField Foam::basicGeochemicalModel::rhol() const
{

    volScalarField rhol_(0.0*Y_[0]);
    forAll(Y_,s)
    {
        rhol_ = rhol_ + Y_[s];
    }

    return rhol_;
}
*/
// ************************************************************************* //
