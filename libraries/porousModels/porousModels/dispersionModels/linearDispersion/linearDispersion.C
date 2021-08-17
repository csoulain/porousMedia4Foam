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

#include "linearDispersion.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dispersionModels
{
    defineTypeNameAndDebug(linearDispersion, 0);

    addToRunTimeSelectionTable
    (
        dispersionModel,
        linearDispersion,
        dictionary
    );
}
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dispersionModels::linearDispersion::linearDispersion
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    dispersionModel(mesh, dict),
    linearDispersionDict_(dict.subDict(typeName+"Coeffs")),
    epsName_(linearDispersionDict_.lookupOrDefault<word>("eps", "eps")),
    UName_(linearDispersionDict_.lookupOrDefault<word>("U", "U")),
    Di_( linearDispersionDict_.lookup("Di") ),
    alphaL_( linearDispersionDict_.lookup("alphaL") ),
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
      dimensionedScalar("Deff",dimensionSet(0, 2, -1, 0, 0),0.0),
      "zeroGradient"
    ),
    eps_(mesh.lookupObject<volScalarField>(epsName_)),
    U_(mesh.lookupObject<volVectorField>(UName_)),
    n_( linearDispersionDict_.lookupOrDefault("n",0.))
{

}

// * * * * * * * * * * * * * * member functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::dispersionModels::linearDispersion::effectiveDispersion() const
{
      return Deff_;
}

void Foam::dispersionModels::linearDispersion::updateDispersion()
{
      Deff_= Foam::pow(eps_,n_)*Di_*(1.+alphaL_/Di_*mag(U_));
}

// alphaL should be zero if eps=1.

// -------------------------------------------------------------------------//
