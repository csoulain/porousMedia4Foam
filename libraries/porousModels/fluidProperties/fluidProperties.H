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
    Foam::fluidProperties

Description
    fluidProperties module for fluid properties in pM4F including
      - density models
      - viscosity models
      - ... models


SourceFiles
    fluidProperties.C

Authors
    - 12/9/2022 C. Soulaine

\*---------------------------------------------------------------------------*/

#ifndef fluidProperties_H
#define fluidProperties_H

#include "volFields.H"
#include "surfaceInterpolate.H"
#include "densityModel.H"
//#include "autoPtr.H"
//#include "viscosity.H"
#include "viscosityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class twoPhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class fluidProperties
{


  protected:

    // Protected data

    const fvMesh &mesh_;

    word fluidName_;

    dictionary fluidDict_;

    autoPtr<densityModel> densityModelPtr_;

    autoPtr<viscosityModel> viscosityModelPtr_;

//    autoPtr<surfaceAreaModel> surfaceAreaModelPtr_;


  // Private Member Functions

    //- Disallow copy construct
    fluidProperties(const fluidProperties&);


  public:

  // Constructors

    //- Construct from components

    fluidProperties
    (
        const fvMesh& mesh,
        const dictionary& dict
    );
/*
    fluidProperties
    (
        const fvMesh& mesh,
        const word & name,
        const dictionary& dict
    );
*/
    //- Destructor
    ~fluidProperties()
    {}


    // Member Functions

    //- Return the composition of the multi-component mixture
//    virtual basicMultiComponentMixture& composition() = 0;

    //- Return the composition of the multi-component mixture
//    virtual const basicMultiComponentMixture& composition() const = 0;

    tmp<volScalarField> rho() const
    {
        return densityModelPtr_->rho();
    }

    void updateDensity()
    {
          densityModelPtr_->updateDensity();
//        return absolutePermeabilityModelPtr_->updatePermeability();
    }

    tmp<volScalarField> nu() const
    {
        return viscosityModelPtr_->nu();
    }

    void updateViscosity()
    {
          viscosityModelPtr_->updateViscosity();
//        return absolutePermeabilityModelPtr_->updatePermeability();
    }

    tmp<volScalarField> mu() const
    {
        return (this->nu())*(this->rho());
    }

/*
    const viscosityModel & viscosity()
    {
        return viscosityModelPtr_();
    }
*/
/*
    viscosity viscosity() const
    {
        return viscosityModelPtr_->nu();
    }
*/

    void update()
    {
        updateDensity();
        updateViscosity();
//        updateSurfaceArea();
    }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
