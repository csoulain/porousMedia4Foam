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

SourceFiles
    porousModel.C

Authors
    - 4/2/2020 C. Soulaine

\*---------------------------------------------------------------------------*/

#ifndef porousModel_H
#define porousModel_H

#include "volFields.H"
#include "surfaceInterpolate.H"
//#include "speciesTable.H"
//#include "runTimeSelectionTables.H"
//#include "surfaceAreaModel.H"
#include "absolutePermeabilityModel.H"
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

    dictionary porousMediaDict_;

    volScalarField eps_;

//    autoPtr<surfaceAreaModel> surfaceAreaModelPtr_;

//    tmp<volScalarField> Ae_;
//    volScalarField Ae_;

    autoPtr<absolutePermeabilityModel> absolutePermeabilityModelPtr_;
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
//        const volVectorField& U,
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

/*

    virtual void updatePorosity();

    virtual void updatePermeability();
*/

/*
    tmp<volScalarField> surfaceArea() const
//    volScalarField surfaceArea() const
    {
        return surfaceAreaModelPtr_->surfaceArea();
    }

    void updateSurfaceArea()
    {
        return surfaceAreaModelPtr_->updateSurfaceArea();
    }
*/

    tmp<volScalarField> absolutePermeability() const
    {
        return absolutePermeabilityModelPtr_->absolutePermeability();
    }

    tmp<volScalarField> invK() const
    {
        dimensionedScalar smallK ("smallK",dimensionSet(0,2,0,0,0,0,0),SMALL);
        return 1./(absolutePermeabilityModelPtr_->absolutePermeability() + smallK);
    }

    surfaceScalarField Kf() const
    {
        surfaceScalarField Kf_
        (
          "Kf",
          fvc::interpolate
          (
              absolutePermeabilityModelPtr_->absolutePermeability(),
              "K"
          )
        );

        return Kf_;
      }
//    tmp<volScalarField> inversePermeability() const
//    {
//        return permeabilityModelPtr_->inversePermeability();
//    }

//    void updatePermeability()
//    {
//        return permeabilityModelPtr_->updatePermeability();
//    }

    const volScalarField & porosity() const
    {
        return eps_;
    }

    void update()
    {}

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
