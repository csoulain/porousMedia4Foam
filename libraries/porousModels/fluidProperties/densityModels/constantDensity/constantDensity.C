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

#include "constantDensity.H"
#include "addToRunTimeSelectionTable.H"


//#include "fvCFD.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace densityModels
{
    defineTypeNameAndDebug(constantDensity, 0);

    addToRunTimeSelectionTable
    (
        densityModel,
        constantDensity,
        dictionary
    );
}
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::densityModels::constantDensity::constantDensity
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    densityModel(mesh, dict),
    constantDensityDict_(dict.subDict(typeName+"Coeffs")),
    rho0_(constantDensityDict_.lookup("rho0"))/*,
    rho_
    (
        IOobject
        (
          "rho",
          mesh_.time().timeName(),
          mesh_,
          IOobject::READ_IF_PRESENT,
          IOobject::NO_WRITE
        ),
        mesh_,
        rho0_,
        "zeroGradient"
    )*/
{
    rho_ = oneField()*rho0_;
//    rho_.correctBoundaryConditions();
}

// * * * * * * * * * * * * * * member functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::densityModels::constantDensity::rho() const
{
//    Info << "fonction rho() de constantDensity1" << nl <<endl;
      return rho_;
}

void Foam::densityModels::constantDensity::updateDensity()
{
    //do nothing

}

// -------------------------------------------------------------------------//
