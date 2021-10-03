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
/*
    updateFromInitialValue_
    (
        KozenyCarmanDict_.lookupOrDefault<Switch>("updateFromInitialValue",false)
    ),*/
    modifiedKozenyCarman_
    (
	KozenyCarmanDict_.lookupOrDefault<Switch>("modifiedKozenyCarman",false)
    ),
    updateFromInitialPoroPerm_
    (
        KozenyCarmanDict_.lookupOrDefault<Switch>("updateFromInitialPoroPerm",false)
    ),
    K_
    (
        IOobject
        (
          "K",
          mesh_.time().timeName(),
          mesh_,
          IOobject::READ_IF_PRESENT,
          IOobject::AUTO_WRITE
        ),
        mesh_,
        K0_,
        //calculatedFvPatchScalarField::typeName
        "zeroGradient"
    ),
    initK_
    (
        IOobject
        (
            "initK",
            mesh_.time().timeName(),
            mesh,
            IOobject::NO_READ
        ),
        K_
    ),
    extrapolateKOnPatchn_( KozenyCarmanDict_.lookupOrDefault("n", 0.66)),
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
    eps_(mesh.lookupObject<volScalarField>(epsName_)),
    eps0_
    (
        IOobject
        (
          "eps0",
          mesh_.time().timeName(),
          mesh_,
          IOobject::READ_IF_PRESENT,
          IOobject::NO_WRITE
        ),
        mesh_,
        1.,
        "zeroGradient"
    )
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
	if (mesh_.time().timeIndex() == 1)
	{
  	 eps0_ = eps_.oldTime();
	}

  if(modifiedKozenyCarman_)
  {
     Info << "Using modified KC" << nl;
     invK_ = Foam::pow(eps0_/(eps_+SMALL),3)/initK_;
  }
  else
  {
   if(updateFromInitialPoroPerm_)
   {
//     Info << "Using Jena KC" << nl;
     //invK_ = Foam::pow(eps0_/(eps_+SMALL),3)/initK_*Foam::pow((1-eps_)/(1-eps0_+SMALL),2);
     invK_ = Foam::pow((1.-eps_),2)/(Foam::pow(eps_+SMALL,3))/K0_;
   }
   else
   {
//    Info << "Ming KC" << nl;
    invK_= Foam::pow((1.-eps_),2)/(Foam::pow(eps_+SMALL,3))
          /(Foam::pow((1.-eps_.oldTime()+SMALL),2)/(Foam::pow(eps_.oldTime()+SMALL,3)))
          *invK_.oldTime();
   }
  }
/*
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
*/
  invK_.max(0.0);



  dimensionedScalar smallInvK_("smallInvK",dimensionSet(0,-2,0,0,0),SMALL);

  K_ = 1./(invK_+smallInvK_);
  /*Extrapolating K_ value on the boundary faces*/

/*
  label patchInlet = mesh_.boundaryMesh().findPatchID("inlet");
  if (patchInlet == -1)
  {
   FatalError << "Inlet patch not found." << exit(FatalError);
  }
*/
  /* SaPa comment - outlet patch not required for extrapolation */
/*
  label patchOutlet = mesh_.boundaryMesh().findPatchID("outlet");
  if (patchOutlet == -1)
  {
   FatalError << "Outlet patch not found." << exit(FatalError);
  }
*/
  /*SaPa comment ends here*/

/*
  label patchWalls = mesh_.boundaryMesh().findPatchID("walls");
  if (patchWalls == -1)
  {
   FatalError << "Walls patch not found." << exit(FatalError);
  }
*/

/*  forAll(mesh_.boundary()[patchInlet],faceI)
  {
   if (eps_[mesh_.boundary()[patchInlet].faceCells()[faceI]] - eps_.oldTime()[mesh_.boundary()[patchInlet].faceCells()[faceI]] >= 0.)
    {
     K_.boundaryFieldRef()[patchInlet][faceI] =
     K_[mesh_.boundary()[patchInlet].faceCells()[faceI]] +
    (Foam::pow(Foam::mag(eps_[mesh_.boundary()[patchInlet].faceCells()[faceI]] - eps_.oldTime()[mesh_.boundary()[patchInlet].faceCells()[faceI]]),extrapolateKOnPatchn_)*K_[mesh_.boundary()[patchInlet].faceCells()[faceI]]);
    }

    else
    {
     K_.boundaryFieldRef()[patchInlet][faceI] =
     K_[mesh_.boundary()[patchInlet].faceCells()[faceI]] -
    (Foam::pow(Foam::mag(eps_[mesh_.boundary()[patchInlet].faceCells()[faceI]]-eps_.oldTime()[mesh_.boundary()[patchInlet].faceCells()[faceI]]),extrapolateKOnPatchn_)*K_[mesh_.boundary()[patchInlet].faceCells()[faceI]]);
    }
  }

  forAll(mesh_.boundary()[patchOutlet],faceI)
  {
   if (eps_[mesh_.boundary()[patchOutlet].faceCells()[faceI]] - eps_.oldTime()[mesh_.boundary()[patchOutlet].faceCells()[faceI]] >= 0.)
   {
    K_.boundaryFieldRef()[patchOutlet][faceI] =
    K_[mesh_.boundary()[patchOutlet].faceCells()[faceI]] +
    (Foam::pow(Foam::mag(eps_[mesh_.boundary()[patchOutlet].faceCells()[faceI]]-eps_.oldTime()[mesh_.boundary()[patchOutlet].faceCells()[faceI]]),extrapolateKOnPatchn_)*K_[mesh_.boundary()[patchOutlet].faceCells()[faceI]]);
   }

   else
   {
    K_.boundaryFieldRef()[patchOutlet][faceI] =
    K_[mesh_.boundary()[patchOutlet].faceCells()[faceI]] -
    (Foam::pow(Foam::mag(eps_[mesh_.boundary()[patchOutlet].faceCells()[faceI]]-eps_.oldTime()[mesh_.boundary()[patchOutlet].faceCells()[faceI]]),extrapolateKOnPatchn_)*K_[mesh_.boundary()[patchOutlet].faceCells()[faceI]]);
   }
  }
*/
/*
  forAll(mesh_.boundary()[patchWalls],faceI)
  {
   if (eps_[mesh_.boundary()[patchWalls].faceCells()[faceI]] - eps_.oldTime()[mesh_.boundary()[patchWalls].faceCells()[faceI]] >= 0.)
   {
    K_.boundaryFieldRef()[patchWalls][faceI] =
    K_[mesh_.boundary()[patchWalls].faceCells()[faceI]] +
    (Foam::pow(Foam::mag(eps_[mesh_.boundary()[patchWalls].faceCells()[faceI]]-eps_.oldTime()[mesh_.boundary()[patchWalls].faceCells()[faceI]]),extrapolateKOnPatchn_)*K_[mesh_.boundary()[patchWalls].faceCells()[faceI]]);
   }

   else
   {
    K_.boundaryFieldRef()[patchWalls][faceI] =
    K_[mesh_.boundary()[patchWalls].faceCells()[faceI]] -
    (Foam::pow(Foam::mag(eps_[mesh_.boundary()[patchWalls].faceCells()[faceI]]-eps_.oldTime()[mesh_.boundary()[patchWalls].faceCells()[faceI]]),extrapolateKOnPatchn_)*K_[mesh_.boundary()[patchWalls].faceCells()[faceI]]);
   }
  }
*/

  Kf_ = fvc::interpolate(K_,"harmonic");

}
// -------------------------------------------------------------------------//
