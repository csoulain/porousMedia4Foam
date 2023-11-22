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

#include "testDensity.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace geochemicalModels
    {
        defineTypeNameAndDebug(testDensity, 0);

        addToRunTimeSelectionTable
        (
            basicGeochemicalModel,
            testDensity,
            dictionary
        );
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::geochemicalModels::testDensity::testDensity
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
      basicGeochemicalModel(mesh, dict),
      densitymodelType_(dict.subDict("fluidProperties").lookup("densityModel")),
      rho_
      (
          IOobject
          (
              "rho",
              mesh.time().timeName(),
              mesh,
              IOobject::READ_IF_PRESENT,
              IOobject::AUTO_WRITE
          ),
          mesh,
          dimensionedScalar("rho",dimensionSet(1,-3,0,0,0,0,0),1e3),
          "zeroGradient"
      )
{
    Y_.resize(0);
}



// * * * * * * * * * * * * * * * * Functions      * * * * * * * * * * * * * * //
Foam::tmp<Foam::volScalarField>
Foam::geochemicalModels::testDensity::rho() const
{
      if(densitymodelType_ == "fromPhreeqc")
      {
            return rho_;
      }
      else
      {
            return Foam::basicGeochemicalModel::rho();
      }
}



void Foam::geochemicalModels::testDensity::updateDensity()
{

//    volScalarField &rho = this->rho();

    if(densitymodelType_ == "fromPhreeqc")
    {
        Info << "Je suis ici " << nl << endl;
        forAll(rho_,cellI)
        {
            rho_[cellI] = cellI;
        }
    }
    else
    {
          Foam::basicGeochemicalModel::updateDensity();
    }

    //densityModelPtr_->updateDensity();

}


// -------------------------------------------------------------------------//


/*
Foam::volScalarField Foam::testDensity::dMl() const
{

    volScalarField dMl_(0.0*fvc::ddt(Y_[0])/this->rhol());
    forAll(Y_,s)
    {
        dMl_ = dMl_ + fvc::ddt(Y_[s])/this->rhol();
    }

    return dMl_;
}
*/

// ************************************************************************* //
