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


    tmp<volScalarField> Se() const
    {
        return basicUnsaturatedGeochemicalModelPtr_->Se();
    }

    tmp<volScalarField> kra() const
    {
        return basicUnsaturatedGeochemicalModelPtr_->kra();
    }

    tmp<volScalarField> krb() const
    {
        return basicUnsaturatedGeochemicalModelPtr_->krb();
    }

    tmp<surfaceScalarField> kraf() const
    {
        return basicUnsaturatedGeochemicalModelPtr_->kraf();
    }

    tmp<surfaceScalarField> krbf() const
    {
        return basicUnsaturatedGeochemicalModelPtr_->krbf();

    }

    tmp<surfaceScalarField> Maf() const
    {
        return basicUnsaturatedGeochemicalModelPtr_->Maf() ;

    }

    tmp<surfaceScalarField> Mbf() const
    {
        return basicUnsaturatedGeochemicalModelPtr_->Mbf() ;

    }

    tmp<surfaceScalarField> Mf() const
    {
        return basicUnsaturatedGeochemicalModelPtr_->Mf() ;

    }

    tmp<surfaceScalarField> Laf() const
    {
        return basicUnsaturatedGeochemicalModelPtr_->Laf() ;

    }

    tmp<surfaceScalarField> Lbf() const
    {
        return basicUnsaturatedGeochemicalModelPtr_->Lbf() ;

    }

    tmp<surfaceScalarField> Lf() const
    {
        return basicUnsaturatedGeochemicalModelPtr_->Lf() ;

    }

    tmp<volScalarField> dkradS() const
    {
        return basicUnsaturatedGeochemicalModelPtr_->dkradS();
    }

    tmp<volScalarField> dkrbdS() const
    {
        return basicUnsaturatedGeochemicalModelPtr_->dkrbdS();
    }

    tmp<volScalarField> pc() const
    {
        return basicUnsaturatedGeochemicalModelPtr_->pc();
    }

    tmp<volScalarField> dpcdS() const
    {
        return basicUnsaturatedGeochemicalModelPtr_->dpcdS();
    }

    Switch activateCapillarity() const
    {
        return basicUnsaturatedGeochemicalModelPtr_->activateCapillarity();
    }


    void update()
    {
        basicUnsaturatedGeochemicalModelPtr_->updateFluidComposition();
        basicUnsaturatedGeochemicalModelPtr_->updateMineralDistribution();
        basicUnsaturatedGeochemicalModelPtr_->updatePorousAssemblageProperties();
        basicUnsaturatedGeochemicalModelPtr_->updatedMinvdRho();
        basicUnsaturatedGeochemicalModelPtr_->updateMobilities();  /// attention on n'update pas pc pour le moment
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
