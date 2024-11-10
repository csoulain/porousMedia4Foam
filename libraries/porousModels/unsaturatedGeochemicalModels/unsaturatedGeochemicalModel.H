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
    Foam::unsaturatedGeochemicalModel

Description
    unsaturatedGeochemicalModel module for geometric porous media properties including
      - porosity
      - absolute permeability models
      - surface area models


SourceFiles
    unsaturatedGeochemicalModel.C

Authors
    - 18/2/2020 C. Soulaine

\*---------------------------------------------------------------------------*/

#ifndef unsaturatedGeochemicalModel_H
#define unsaturatedGeochemicalModel_H

#include "volFields.H"
#include "surfaceInterpolate.H"
//#include "speciesTable.H"
//#include "runTimeSelectionTables.H"
//#include "fluidProperties.H"
#include "basicUnsaturatedGeochemicalModel.H"
//#include "autoPtr.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class twoPhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class unsaturatedGeochemicalModel
{


  protected:

    // Protected data

    const fvMesh &mesh_;

//    dictionary unsaturatedGeochemicalModelDict_;

//    dictionary fluidPropertiesDict_;


    autoPtr<basicUnsaturatedGeochemicalModel> basicUnsaturatedGeochemicalModelPtr_;


  // Private Member Functions

    //- Disallow copy construct
    unsaturatedGeochemicalModel(const unsaturatedGeochemicalModel&);


  public:

  // Constructors

    //- Construct from components

    unsaturatedGeochemicalModel
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const volScalarField &Sb,
        const incompressiblePhase &phasea,
        const incompressiblePhase &phaseb
    );


    //- Destructor
    ~unsaturatedGeochemicalModel()
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
        return basicUnsaturatedGeochemicalModelPtr_->absolutePermeability();
    }

    tmp<volScalarField> inversePermeability() const
    {
        return basicUnsaturatedGeochemicalModelPtr_->inversePermeability();
    }

    const surfaceScalarField & Kf() const
    {
        return basicUnsaturatedGeochemicalModelPtr_->Kf();
    }

    const volScalarField & porosity()  const
    {
        return basicUnsaturatedGeochemicalModelPtr_->porosity();
    }


    const volScalarField & dMinvdRho()  const
    {
        return basicUnsaturatedGeochemicalModelPtr_->dMinvdRho();
    }

    const volTensorField & effectiveDispersionTensor()  const
    {
        return basicUnsaturatedGeochemicalModelPtr_->effectiveDispersionTensor();
    }


    const volScalarField & rho()  const
    {
        return basicUnsaturatedGeochemicalModelPtr_->rho();
    }

    const surfaceScalarField & rhof()  const
    {
        return fvc::interpolate(basicUnsaturatedGeochemicalModelPtr_->rho());
    }

/*
    const viscosityModel & viscosity()  const
    {
        return basicUnsaturatedGeochemicalModelPtr_->viscosity();
    }
*/
    const volScalarField & nu()  const
    {
        return basicUnsaturatedGeochemicalModelPtr_->nu();
    }

    tmp<volScalarField> mu() const
    {

        dimensionedScalar nu00 ("nu", dimensionSet(1,-1,-1,0,0,0,0), 0.0);

        volScalarField mu
        (
            IOobject
            (
              "nu",
              mesh_.time().timeName(),
              mesh_,
              IOobject::READ_IF_PRESENT,
              IOobject::NO_WRITE
            ),
            mesh_,
            nu00,
            "zeroGradient"
        );
        mu = (this->nu())*(this->rho());
        mu.correctBoundaryConditions();


        return mu;
    }
/*
    const viscosityModel & viscosity()
    {
        return basicUnsaturatedGeochemicalModelPtr_->viscosity();
    }
*/
    void update()
    {
        basicUnsaturatedGeochemicalModelPtr_->updateFluidComposition();
        basicUnsaturatedGeochemicalModelPtr_->updateMineralDistribution();
        basicUnsaturatedGeochemicalModelPtr_->updatePorosity();
        basicUnsaturatedGeochemicalModelPtr_->updateSurfaceArea();
        basicUnsaturatedGeochemicalModelPtr_->updatePermeability();
        basicUnsaturatedGeochemicalModelPtr_->updateDispersionTensor();
        basicUnsaturatedGeochemicalModelPtr_->updatedMinvdRho();
        basicUnsaturatedGeochemicalModelPtr_->updateDensity();
        basicUnsaturatedGeochemicalModelPtr_->updateViscosity();
        Info << "geochemistry updated" <<nl <<endl;
    }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
