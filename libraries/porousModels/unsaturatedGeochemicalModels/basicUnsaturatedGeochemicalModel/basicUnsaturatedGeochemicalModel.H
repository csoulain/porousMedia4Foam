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
    Foam::basicUnsaturatedGeochemicalModel

Description
    basicUnsaturatedGeochemicalModel module

SourceFiles
    basicUnsaturatedGeochemicalModel.C

Authors
    17/02/2020 : C. Soulaine


\*---------------------------------------------------------------------------*/

#ifndef basicUnsaturatedGeochemicalModel_H
#define basicUnsaturatedGeochemicalModel_H

//#include "dictionary.H"
//#include "volFieldsFwd.H"
#include "volFields.H"
#include "speciesTable.H"
#include "runTimeSelectionTables.H"
#include "unsaturatedPorousModel.H"
#include "fluidProperties.H"
//#include "autoPtr.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class twoPhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class basicUnsaturatedGeochemicalModel : public fluidProperties
{


  protected:

    // Protected data

    const fvMesh &mesh_;

    dictionary geochemicalModelDict_;

    dictionary fluidPropertiesDict_;

    //fluidProperties fluidProperties_;


    speciesTable mineralList_;

    PtrList<volScalarField> Ys_;

    volScalarField inertMineral_;

    volScalarField eps_;
    volScalarField eps0_;

    PtrList<dimensionedScalar> rhos_;


    volScalarField dMinvdRho_;   /// mass transfer from solid to liquid

    volScalarField dGasvdRho_;   /// mass transfer from gas to liquid

    PtrList<porousModel> porousMedia_;


//    autoPtr<densityModel> densityModelPtr_;

//    autoPtr<viscosityModel> viscosityModelPtr_;

    autoPtr<absolutePermeabilityModel> absolutePermeabilityModelPtr_;

    autoPtr<dispersionTensorModel> dispersionTensorModelPtr_;


    autoPtr<reducedSaturationModel> reducedSaturationModelPtr_;

    autoPtr<relativePermeabilityModel> relativePermeabilityModelPtr_;

    autoPtr<capillarityModel> capillarityModelPtr_;

    autoPtr<incompressiblePhase> phasea;
    autoPtr<incompressiblePhase> phaseb;
    dimensionedScalar rhoa_;
    dimensionedScalar rhob_;	
    //ORIG->
    //word phiName_;
    word phiaName_;
    word phibName_;
    const surfaceScalarField &phia_;
    const surfaceScalarField &phib_;


    word phiName_;

    const surfaceScalarField &phi_;

    word SbName_;
 
    const volScalarField &Sb_;


/*
    // concentration fields
    PtrList<volScalarField> Y_;

    PtrList<volScalarField> Ya_; //allocate components in phase a
    PtrList<volScalarField> Yb_; //allocate components in phase b
*/

  // Private Member Functions

    //- Disallow copy construct
    basicUnsaturatedGeochemicalModel(const basicUnsaturatedGeochemicalModel&);


  public:

    //- Runtime type information
      TypeName("basicUnsaturatedGeochemicalModel");


    // Declare run-time constructor selection table

      declareRunTimeSelectionTable
      (
          autoPtr,
          basicUnsaturatedGeochemicalModel,
          dictionary,
          (
              const fvMesh& mesh,
              const dictionary& dict
          ),
          (mesh, dict)
      );

    // Selectors

        //- Return a reference to the selected viscosity model
        static autoPtr<basicUnsaturatedGeochemicalModel> New
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


  // Constructors

    //- Construct from components

    basicUnsaturatedGeochemicalModel
    (
        const fvMesh& mesh,
        const dictionary& dict
    );

    //- Destructor
    virtual ~basicUnsaturatedGeochemicalModel()
    {}


    // Member Functions

    virtual void updatePorosity();

  //  virtual volScalarField dMs() const = 0;

    tmp<volScalarField> porosity() const
    {
        return eps_;
    }

    tmp<volScalarField> absolutePermeability() const
    {
        return absolutePermeabilityModelPtr_->absolutePermeability();
    }

    tmp<volScalarField> inversePermeability() const
    {
        return absolutePermeabilityModelPtr_->inversePermeability();
    }

    tmp<surfaceScalarField> Kf() const
    {
        return absolutePermeabilityModelPtr_->Kf();
    }

    void updatePermeability()
    {
        absolutePermeabilityModelPtr_->updatePermeability();
    }

    void updateSurfaceArea()
    {
        forAll(porousMedia_,s)
        {
            porousMedia_[s].updateSurfaceArea();
        }
    }

    tmp<volTensorField> effectiveDispersionTensor() const
    {
        return dispersionTensorModelPtr_->effectiveDispersionTensor();
    }

    void updateDispersionTensor()
    {
        dispersionTensorModelPtr_->updateDispersionTensor();
    }


    tmp<volScalarField> kra() const
    {
        return relativePermeabilityModelPtr_->kra();
    }
    tmp<surfaceScalarField> kraf() const
    {
    	return relativePermeabilityModelPtr_->kraf();
    }

    tmp<volScalarField> krb() const
    {
        return relativePermeabilityModelPtr_->krb();
    }
    tmp<surfaceScalarField> krbf() const
    {
	return relativePermeabilityModelPtr_->krbf();
    }

    tmp<volScalarField> dkradS() const
    {
    	return relativePermeabilityModelPtr_->dkradS();
    }
    tmp<volScalarField> dkrbdS() const
    {
     	return relativePermeabilityModelPtr_->dkrbdS();
    }
    /*
    Switch activateCapillarity() const
    {
    	return capillarityModelPtr_->activateCapillarity();
    }
    */
    tmp<volScalarField> dpcdS() const
    {
     	return capillarityModelPtr_->dpcdS();
    }

/*
    virtual tmp<volScalarField> rho() const
    {
        return fluidProperties.rho();
    }

    virtual void updateDensity()
    {
        this->updateDensity();

    }
    */
/*

    virtual tmp<volScalarField> nu() const
    {
//        return viscosityModelPtr_->nu();
        return fluidProperties_.nu();
    }

    virtual void updateViscosity()
    {
//        viscosityModelPtr_->updateViscosity();
        fluidProperties_.updateViscosity();
    }
*/
/*
    const viscosityModel & viscosity()
    {
        return viscosityModelPtr_();
    }
*/

    tmp<volScalarField> dMinvdRho() const /// solid to liquid
    {
        return dMinvdRho_;
    }

    virtual void updatedMinvdRho();


    tmp<volScalarField> dGasvdRho() const  // gas to liquid
    {
     return dGasvdRho_;
    }

    virtual void updatedGasvdRho();


    virtual void updateFluidComposition() = 0;

    virtual void updateMineralDistribution() = 0;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
