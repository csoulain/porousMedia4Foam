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

#include "taylorAris.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dispersionModels
{
    defineTypeNameAndDebug(taylorAris, 0);

    addToRunTimeSelectionTable
    (
        dispersionModel,
        taylorAris,
        dictionary
    );
}
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dispersionModels::taylorAris::taylorAris
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    dispersionModel(mesh, dict),
    taylorArisDict_(dict.subDict(typeName+"Coeffs")),
    UName_(taylorArisDict_.lookupOrDefault<word>("U", "U")),
    Di_( taylorArisDict_.lookup("Di") ),
    Deff_
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
    U_(mesh.lookupObject<volVectorField>(UName_)),
    d_( taylorArisDict_.lookup("channelDiameter") )

{

}

// * * * * * * * * * * * * * * member functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::dispersionModels::taylorAris::effectiveDispersion() const
{
      return Deff_;
}

void Foam::dispersionModels::taylorAris::updateDispersion()
{
      Deff_= Di_*(1 + Foam::pow(d_*0.5,2)*Foam::pow(mag(U_),2)/(48*Foam::pow(Di_,2)));
}
// -------------------------------------------------------------------------//