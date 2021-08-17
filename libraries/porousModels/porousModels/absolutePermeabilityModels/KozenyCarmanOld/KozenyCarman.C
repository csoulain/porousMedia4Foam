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

#include "KozenyCarman.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace absolutePermeabilityModels
{
    defineTypeNameAndDebug(KozenyCarman, 0);

    addToRunTimeSelectionTable
    (
        absolutePermeabilityModel,
        KozenyCarman,
        dictionary
    );
}
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::absolutePermeabilityModels::KozenyCarman::KozenyCarman
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    absolutePermeabilityModel(mesh, dict),
    KozenyCarmanDict_(dict.subDict(typeName+"Coeffs")),
    epsName_(KozenyCarmanDict_.lookupOrDefault<word>("eps", "eps")),
    K0_(KozenyCarmanDict_.lookupOrDefault
    (
        "K0",
        dimensionedScalar("K0",dimensionSet(0,2,0,0,0,0,0),SMALL))
    ),
    updateFromInitialValue_
    (
        KozenyCarmanDict_.lookupOrDefault<Switch>("updateFromInitialValue",false)
    ),
    K_
    (
        IOobject
        (
          "K",
          mesh_.time().timeName(),
          mesh_,
          IOobject::READ_IF_PRESENT,
          IOobject::NO_WRITE
        ),
        mesh_,
        K0_,
        "zeroGradient"
    ),
    invK_
    (
        IOobject
        (
          "invK",
          mesh_.time().timeName(),
          mesh_,
          IOobject::NO_READ,
          IOobject::NO_WRITE
        ),
        1/K_,
        "zeroGradient"
    ),
    Kf_("Kf", fvc::interpolate(K_,"K")),
    eps_(mesh.lookupObject<volScalarField>(epsName_))
{

//    updatePermeability();

}

// * * * * * * * * * * * * * * member functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::absolutePermeabilityModels::KozenyCarman::absolutePermeability() const
{
      return K_;
}

Foam::tmp<Foam::volScalarField>
Foam::absolutePermeabilityModels::KozenyCarman::inversePermeability() const
{
      return invK_;
}

Foam::tmp<Foam::surfaceScalarField>
Foam::absolutePermeabilityModels::KozenyCarman::Kf() const
{
      return Kf_;
}


void Foam::absolutePermeabilityModels::KozenyCarman::updatePermeability()
{

  if(updateFromInitialValue_)
  {
      invK_= Foam::pow((1.-eps_),2)/(Foam::pow(eps_+SMALL,3))/K0_;
  }
  else
  {
    invK_= Foam::pow((1.-eps_),2)/(Foam::pow(eps_+SMALL,3))
          /(Foam::pow((1.-eps_.oldTime()-SMALL),2)/(Foam::pow(eps_.oldTime()+SMALL,3)))
          *invK_.oldTime();
  }

  invK_.max(0.0);



  dimensionedScalar smallInvK_("smallInvK",dimensionSet(0,-2,0,0,0),SMALL);

  K_ = 1./(invK_+smallInvK_);

  Kf_ = fvc::interpolate(K_,"K");

}
// -------------------------------------------------------------------------//
