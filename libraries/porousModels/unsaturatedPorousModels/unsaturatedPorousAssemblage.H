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

Class
    Foam::unsaturatedPorousAssemblage

Description
    unsaturatedPorousAssemblage module for geometric porous media properties including
      - porosity
      - absolute permeability models
      - surface area models


SourceFiles
    unsaturatedPorousAssemblage.C

Authors
    - 4/2/2020 C. Soulaine

\*---------------------------------------------------------------------------*/

#ifndef unsaturatedPorousAssemblage_H
#define unsaturatedPorousAssemblage_H

#include "volFields.H"
#include "surfaceInterpolate.H"
#include "porousAssemblage.H"
#include "incompressiblePhase.H"
#include "reducedSaturationModel.H"
#include "relativePermeabilityModel.H"
#include "capillarityModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class twoPhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class unsaturatedPorousAssemblage : public porousAssemblage
{


  protected:

    // Protected data

    const fvMesh &mesh_;

    dictionary porousMediaDict_;

    const incompressiblePhase &phasea_;

    const incompressiblePhase &phaseb_;

    const volScalarField &Sb_;

    autoPtr<reducedSaturationModel> reducedSaturationModelPtr_;

    autoPtr<relativePermeabilityModel> relativePermeabilityModelPtr_;

    autoPtr<capillarityModel> capillarityModelPtr_;

    surfaceScalarField Maf_,Mbf_,Mf_;
    surfaceScalarField Laf_,Lbf_,Lf_;


  // Private Member Functions

    //- Disallow copy construct
    unsaturatedPorousAssemblage(const unsaturatedPorousAssemblage&);


  public:

  // Constructors

    //- Construct from components

    unsaturatedPorousAssemblage
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const volScalarField &Sb,
        const incompressiblePhase &phasea,
        const incompressiblePhase &phaseb
    );

    //- Destructor
    ~unsaturatedPorousAssemblage()
    {}


    // Member Functions

    //- Return the composition of the multi-component mixture
//    virtual basicMultiComponentMixture& composition() = 0;

    //- Return the composition of the multi-component mixture
//    virtual const basicMultiComponentMixture& composition() const = 0;


    // Member Functions
    tmp<volScalarField> Se() const
    {
        return reducedSaturationModelPtr_->Se();
    }

    tmp<volScalarField> kra() const
    {
        return relativePermeabilityModelPtr_->kra();
    }

    tmp<volScalarField> krb() const
    {
        return relativePermeabilityModelPtr_->krb();
    }

    tmp<surfaceScalarField> kraf() const
    {
        return relativePermeabilityModelPtr_->kraf();
    }

    tmp<surfaceScalarField> krbf() const
    {
        return relativePermeabilityModelPtr_->krbf();

    }

    tmp<surfaceScalarField> Maf() const
    {
        return Maf_ ;

    }

    tmp<surfaceScalarField> Mbf() const
    {
        return Mbf_ ;

    }

    tmp<surfaceScalarField> Mf() const
    {
        return Mf_ ;

    }

    tmp<surfaceScalarField> Laf() const
    {
        return Laf_ ;

    }

    tmp<surfaceScalarField> Lbf() const
    {
        return Lbf_ ;

    }

    tmp<surfaceScalarField> Lf() const
    {
        return Lf_ ;

    }

/*
    tmp<surfaceScalarField> Faf() const
    {
        return Faf_ ;

    }

    tmp<surfaceScalarField> Fbf() const
    {
        return Fbf_ ;

    }
*/

    tmp<volScalarField> dkradS() const
    {
        return relativePermeabilityModelPtr_->dkradS();
    }

    tmp<volScalarField> dkrbdS() const
    {
        return relativePermeabilityModelPtr_->dkrbdS();
    }

    tmp<volScalarField> pc() const
    {
        return capillarityModelPtr_->pc();
    }

    tmp<volScalarField> dpcdS() const
    {
        return capillarityModelPtr_->dpcdS();
    }

    Switch activateCapillarity() const
    {
        return capillarityModelPtr_->activateCapillarity();
    }

    void updateMobilities()
    {
        reducedSaturationModelPtr_->updateReducedSaturation();
        relativePermeabilityModelPtr_->correct();
        Maf_ = this->Kf()*this->kraf()/phasea_.mu();
        Mbf_ = Kf()*krbf()/phaseb_.mu();
        Mf_  = Maf_+Mbf_;
        Laf_ = Maf_*phasea_.rho();
        Lbf_ = Mbf_*phaseb_.rho();
        Lf_  = Laf_+Lbf_;

//        Faf_ = Maf_/Mf_;
//        Fbf_ = Mbf_/Mf_;
    }


/*
    void update()
    {
        updatePermeability();
        updateSurfaceArea();
        updateDispersionTensor();
    }
    */

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
