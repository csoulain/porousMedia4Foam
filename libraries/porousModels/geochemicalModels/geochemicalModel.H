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
    Foam::geochemicalModel

Description
    geochemicalModel module for geometric porous media properties including
      - porosity
      - absolute permeability models
      - surface area models


SourceFiles
    geochemicalModel.C

Authors
    - 18/2/2020 C. Soulaine

\*---------------------------------------------------------------------------*/

#ifndef geochemicalModel_H
#define geochemicalModel_H

#include "volFields.H"
#include "surfaceInterpolate.H"
//#include "speciesTable.H"
//#include "runTimeSelectionTables.H"
//#include "fluidProperties.H"
#include "basicGeochemicalModel.H"
//#include "autoPtr.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class twoPhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class geochemicalModel
{


  protected:

    // Protected data

    const fvMesh &mesh_;

//    dictionary geochemicalModelDict_;

//    dictionary fluidPropertiesDict_;

//    volScalarField Ys_;

//    volScalarField eps_;

    autoPtr<basicGeochemicalModel> basicGeochemicalModelPtr_;

  //  autoPtr<densityModel> densityModelPtr_;

/*
    tmp<volScalarField> invK_;
*/
  // Private Member Functions

    //- Disallow copy construct
    geochemicalModel(const geochemicalModel&);


  public:

  // Constructors

    //- Construct from components

    geochemicalModel
    (
        const fvMesh& mesh,
        const dictionary& dict
    );

    geochemicalModel
    (
        const fvMesh& mesh,
        const volScalarField& Ys,
        const dictionary& dict
    );

    //- Destructor
    ~geochemicalModel()
    {}


    // Member Functions

    //- Return the composition of the multi-component mixture
//    virtual basicMultiComponentMixture& composition() = 0;

    //- Return the composition of the multi-component mixture
//    virtual const basicMultiComponentMixture& composition() const = 0;

/*
    tmp<volScalarField> surfaceArea() const
    {
        return surfaceAreaModelPtr_->surfaceArea();
    }

    void updateSurfaceArea()
    {
        return surfaceAreaModelPtr_->updateSurfaceArea();
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
        return absolutePermeabilityModelPtr_->updatePermeability();
    }

    tmp<surfaceScalarField> Kf() const
    {
        return absolutePermeabilityModelPtr_->Kf();
//          return fvc::interpolate(absolutePermeability(),"K");

    }



*/

    const volScalarField & absolutePermeability()  const
    {
        return basicGeochemicalModelPtr_->absolutePermeability();
    }

    tmp<volScalarField> inversePermeability() const
    {
        return basicGeochemicalModelPtr_->inversePermeability();
    }

    const surfaceScalarField & Kf() const
    {
        return basicGeochemicalModelPtr_->Kf();
    }

    const volScalarField & porosity()  const
    {
        return basicGeochemicalModelPtr_->porosity();
    }


    const volScalarField & dMinvdRho()  const
    {
        return basicGeochemicalModelPtr_->dMinvdRho();
    }

    const volTensorField & effectiveDispersionTensor()  const
    {
        return basicGeochemicalModelPtr_->effectiveDispersionTensor();
    }


    const volScalarField & rho()  const
    {
        return basicGeochemicalModelPtr_->rho();
    }

    const surfaceScalarField & rhof()  const
    {
        return fvc::interpolate(basicGeochemicalModelPtr_->rho());
    }

/*
    const viscosityModel & viscosity()  const
    {
        return basicGeochemicalModelPtr_->viscosity();
    }
*/
    const volScalarField & nu()  const
    {
        return basicGeochemicalModelPtr_->nu();
    }
/*
    const viscosityModel & viscosity()
    {
        return basicGeochemicalModelPtr_->viscosity();
    }
*/
    void update()
    {
        basicGeochemicalModelPtr_->updateFluidComposition();
        basicGeochemicalModelPtr_->updateMineralDistribution();
        basicGeochemicalModelPtr_->updatePorosity();
        basicGeochemicalModelPtr_->updateSurfaceArea();
        basicGeochemicalModelPtr_->updatePermeability();
        basicGeochemicalModelPtr_->updateDispersionTensor();
        basicGeochemicalModelPtr_->updatedMinvdRho();
        basicGeochemicalModelPtr_->updateDensity();
        Info << "geochemistry updated" <<nl <<endl;
    }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
