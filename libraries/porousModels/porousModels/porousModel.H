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
    Foam::porousModel

Description
    porousModel module for geometric porous media properties including
      - porosity
      - absolute permeability models
      - surface area models


SourceFiles
    porousModel.C

Authors
    - 4/2/2020 C. Soulaine

\*---------------------------------------------------------------------------*/

#ifndef porousModel_H
#define porousModel_H

#include "volFields.H"
#include "surfaceInterpolate.H"
#include "absolutePermeabilityModel.H"
//#include "dispersionModel.H"
#include "dispersionTensorModel.H"
#include "surfaceAreaModel.H"
//#include "autoPtr.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class twoPhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class porousModel
{


  protected:

    // Protected data

    const fvMesh &mesh_;

    word porousMediaName_;

    dictionary porousMediaDict_;

    volScalarField Ys_;

    const volScalarField &Yss_; //Cyp temp addition for Ys update
                                //correctly in geochemicalModel

    volScalarField eps_;

    autoPtr<absolutePermeabilityModel> absolutePermeabilityModelPtr_;

//    autoPtr<dispersionModel> dispersionModelPtr_;

    autoPtr<dispersionTensorModel> dispersionTensorModelPtr_;


    autoPtr<surfaceAreaModel> surfaceAreaModelPtr_;


/*
    tmp<volScalarField> invK_;
*/
  // Private Member Functions

    //- Disallow copy construct
    porousModel(const porousModel&);


  public:

  // Constructors

    //- Construct from components

    porousModel
    (
        const fvMesh& mesh,
        const dictionary& dict
    );

    porousModel
    (
        const fvMesh& mesh,
        const word & name,
        const volScalarField& Ys,
        const dictionary& dict
    );

    //- Destructor
    ~porousModel()
    {}


    // Member Functions

    //- Return the composition of the multi-component mixture
//    virtual basicMultiComponentMixture& composition() = 0;

    //- Return the composition of the multi-component mixture
//    virtual const basicMultiComponentMixture& composition() const = 0;


    tmp<volScalarField> surfaceArea() const
    {
        return surfaceAreaModelPtr_->surfaceArea();
    }

    void updateSurfaceArea()
    {
        surfaceAreaModelPtr_->updateSurfaceArea();
//        return surfaceAreaModelPtr_->updateSurfaceArea();
    }


    tmp<volScalarField> absolutePermeability() const
    {
        return absolutePermeabilityModelPtr_->absolutePermeability();
    }

    tmp<volScalarField> invK() const
    {
        return absolutePermeabilityModelPtr_->inversePermeability();
    }

    void updatePermeability()
    {
          absolutePermeabilityModelPtr_->updatePermeability();
//        return absolutePermeabilityModelPtr_->updatePermeability();
    }

    tmp<surfaceScalarField> Kf() const
    {
        return absolutePermeabilityModelPtr_->Kf();
//          return fvc::interpolate(absolutePermeability(),"K");

    }
/*
    tmp<volScalarField> effectiveDispersion() const
    {
        return dispersionModelPtr_->effectiveDispersion();
    }

    void updateDispersion()
    {
        dispersionModelPtr_->updateDispersion();
    }
*/
    tmp<volTensorField> effectiveDispersionTensor() const
    {
        return dispersionTensorModelPtr_->effectiveDispersionTensor();
    }

    void updateDispersionTensor()
    {
        dispersionTensorModelPtr_->updateDispersionTensor();
    }


    const volScalarField & porosity() const
    {
        return eps_;
    }

    void update()
    {
        updatePermeability();
        updateSurfaceArea();
    //    updateDispersion();
        updateDispersionTensor();
    }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
