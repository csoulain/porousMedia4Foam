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

#include "generalDispersionModel.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dispersionModels
{
    defineTypeNameAndDebug(generalDispersionModel, 0);

    addToRunTimeSelectionTable
    (
        dispersionModel,
        generalDispersionModel,
        dictionary
    );
}
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dispersionModels::generalDispersionModel::generalDispersionModel
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    dispersionModel(mesh, dict),
    generalDispersionModelDict_(dict.subDict(typeName+"Coeffs")),
    epsName_(generalDispersionModelDict_.lookupOrDefault<word>("eps", "eps")),
    UName_(generalDispersionModelDict_.lookupOrDefault<word>("U", "U")),
    Di_( generalDispersionModelDict_.lookup("Di") ),
    alphaL_( generalDispersionModelDict_.lookup("alphaL") ),
    Deff
     (
      IOobject
      (
       "Deff",
       mesh_.time().timeName(),
       mesh_,
       IOobject::READ_IF_PRESENT,
       IOobject::AUTO_WRITE
      ),
      mesh_,
      dimensionedScalar("Deff",dimensionSet(0, 2, -1, 0, 0),0.0),
      "zeroGradient"
    ),
    eps_(mesh.lookupObject<volScalarField>(epsName_)),
    U_(mesh.lookupObject<volVectorField>(UName_))
{

}

// * * * * * * * * * * * * * * member functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::dispersionModels::generalDispersionModel::effectiveDispersion() const
{
      return Deff;
}

void Foam::dispersionModels::generalDispersionModel::updateDispersion()
{
      Deff= Foam::pow(eps_,1)*Di_*(1.+alphaL_/Di_*mag(U_));
}
// -------------------------------------------------------------------------//
