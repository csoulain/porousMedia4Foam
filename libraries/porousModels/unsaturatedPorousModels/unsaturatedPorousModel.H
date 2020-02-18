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
    Foam::unsaturatedPorousModel

Description
    unsaturatedPorousModel module for porous media properties
    in unsaturated (two-fluid) conditions including
      - relative permeability models
      - capillary pressure curve models

SourceFiles
    unsaturatedPorousModel.C

Authors
    - 5/2/2020 C. Soulaine

\*---------------------------------------------------------------------------*/

#ifndef unsaturatedunsaturatedPorousModel_H
#define unsaturatedPorousModel_H

#include "porousModel.H"
#include "incompressiblePhase.H"
#include "reducedSaturationModel.H"
#include "relativePermeabilityModel.H"
#include "capillarityModel.H"
//#include "autoPtr.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class twoPhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class unsaturatedPorousModel :  public porousModel
{


  protected:

    // Protected data

    const fvMesh &mesh_;

    dictionary porousMediaDict_;

    const volScalarField &Sb_;

    const incompressiblePhase &phasea_;

    const incompressiblePhase &phaseb_;

    autoPtr<reducedSaturationModel> reducedSaturationModelPtr_;

    autoPtr<relativePermeabilityModel> relativePermeabilityModelPtr_;

    autoPtr<capillarityModel> capillarityModelPtr_;

    surfaceScalarField Maf_,Mbf_,Mf_;
    surfaceScalarField Laf_,Lbf_,Lf_;
//    surfaceScalarField Faf_,Fbf_;


  // Private Member Functions

    //- Disallow copy construct
    unsaturatedPorousModel(const unsaturatedPorousModel&);


  public:

  // Constructors

    //- Construct from components

    unsaturatedPorousModel
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const volScalarField &Sb,
        const incompressiblePhase &phasea,
        const incompressiblePhase &phaseb
    );

    //- Destructor
    ~unsaturatedPorousModel()
    {}


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
        Maf_ = this->Kf()*this->kraf()/phasea_.mu();
        Mbf_ = Kf()*krbf()/phaseb_.mu();
        Mf_  = Maf_+Mbf_;
        Laf_ = Maf_*phasea_.rho();
        Lbf_ = Mbf_*phaseb_.rho();
        Lf_  = Laf_+Lbf_;

//        Faf_ = Maf_/Mf_;
//        Fbf_ = Mbf_/Mf_;
    }

    void update()
    {
        reducedSaturationModelPtr_->updateReducedSaturation();
        relativePermeabilityModelPtr_->correct();
        capillarityModelPtr_->correct();
        this->updateMobilities();
    }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
