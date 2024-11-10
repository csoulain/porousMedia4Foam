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
    Foam::porousAssemblage

Description
    porousAssemblage module for geometric porous media properties including
      - porosity
      - absolute permeability models
      - surface area models


SourceFiles
    porousAssemblage.C

Authors
    - 4/2/2020 C. Soulaine

\*---------------------------------------------------------------------------*/

#ifndef porousAssemblage_H
#define porousAssemblage_H

#include "volFields.H"
#include "surfaceInterpolate.H"
#include "absolutePermeabilityModel.H"
#include "dispersionTensorModel.H"
#include "surfaceAreaModel.H"
//#include "autoPtr.H"
#include "speciesTable.H"
#include "mineralModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class twoPhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class porousAssemblage
{


  protected:

    // Protected data

    const fvMesh &mesh_;

    speciesTable mineralList_;

    PtrList<volScalarField> Ys_;


    volScalarField inertMineral_;

    volScalarField eps_;
    volScalarField eps0_;

    PtrList<dimensionedScalar> rhos_;

    PtrList<volScalarField> mineralSurfaceArea_;


    PtrList<mineralModel> mineral_;


    autoPtr<absolutePermeabilityModel> absolutePermeabilityModelPtr_;

    autoPtr<dispersionTensorModel> dispersionTensorModelPtr_;



  // Private Member Functions

    //- Disallow copy construct
    porousAssemblage(const porousAssemblage&);


  public:

  // Constructors

    //- Construct from components

    porousAssemblage
    (
        const fvMesh& mesh,
        const dictionary& dict
    );

    //- Destructor
    ~porousAssemblage()
    {}


    // Member Functions

    //- Return the composition of the multi-component mixture
//    virtual basicMultiComponentMixture& composition() = 0;

    //- Return the composition of the multi-component mixture
//    virtual const basicMultiComponentMixture& composition() const = 0;


    tmp<volScalarField> absolutePermeability() const
    {
        return absolutePermeabilityModelPtr_->absolutePermeability();
    }

    tmp<volScalarField> inversePermeability() const
    {
        return absolutePermeabilityModelPtr_->inversePermeability();
    }

    void updatePermeability()
    {
          absolutePermeabilityModelPtr_->updatePermeability();
    }

    //tmp<surfaceScalarField> Kf() const
    const surfaceScalarField& Kf() const
    {
        Info << "avant Kf dans porousAssemblage" <<endl;
        return absolutePermeabilityModelPtr_->Kf();
//          return fvc::interpolate(absolutePermeability(),"K");

    }

    tmp<volTensorField> effectiveDispersionTensor() const
    {
        return dispersionTensorModelPtr_->effectiveDispersionTensor();
    }

    void updateDispersionTensor()
    {
        dispersionTensorModelPtr_->updateDispersionTensor();
    }


    void updatePorosity();

    const volScalarField & porosity() const
    {
        return eps_;
    }


    void updateSurfaceArea()
    {
        forAll(mineralList_,s)
        {
            mineral_[s].updateSurfaceArea();
        }
    }


    void update()
    {
        updatePorosity();
        updatePermeability();
        updateSurfaceArea();
        updateDispersionTensor();
    }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
