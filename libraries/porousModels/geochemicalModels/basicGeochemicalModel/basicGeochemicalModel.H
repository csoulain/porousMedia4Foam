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
    Foam::twoPhaseMixture

Description
    basicGeochemicalModel module

SourceFiles
    phreeqc.C

\*---------------------------------------------------------------------------*/

#ifndef basicGeochemicalModel_H
#define basicGeochemicalModel_H

//#include "dictionary.H"
//#include "volFieldsFwd.H"
#include "volFields.H"
#include "speciesTable.H"
#include "runTimeSelectionTables.H"
#include "porousModel.H"
//#include "autoPtr.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class twoPhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class basicGeochemicalModel
{


  protected:

    // Protected data

    const fvMesh &mesh_;

    speciesTable mineralList_;

    PtrList<volScalarField> Ys_;

    volScalarField inertMineral_;

    volScalarField eps_;

    PtrList<porousModel> porousMedia_;

//    tmp<volScalarField> Ae_;
//    volScalarField Ae_;

    autoPtr<absolutePermeabilityModel> absolutePermeabilityModelPtr_;
/*
    tmp<volScalarField> invK_;
*/
  // Private Member Functions

    //- Disallow copy construct
    basicGeochemicalModel(const basicGeochemicalModel&);


  public:

    //- Runtime type information
      TypeName("basicGeochemicalModel");


    // Declare run-time constructor selection table

      declareRunTimeSelectionTable
      (
          autoPtr,
          basicGeochemicalModel,
          dictionary,
          (
              const fvMesh& mesh,
              const dictionary& dict
          ),
          (mesh, dict)
      );

    // Selectors

        //- Return a reference to the selected viscosity model
        static autoPtr<basicGeochemicalModel> New
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


  // Constructors

    //- Construct from components

    basicGeochemicalModel
    (
        const fvMesh& mesh,
        const dictionary& dict
    );

    //- Destructor
    virtual ~basicGeochemicalModel()
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
/*    rm -f $FOAM_USER_LIBBIN/libgeochemicalModels.so

    tmp<volScalarField> surfaceArea() const
//    volScalarField surfaceArea() const
    {
        return surfaceAreaModelPtr_->surfaceArea();
    }

    void updateSurfaceArea()
    {
        return surfaceAreaModelPtr_->updateSurfaceArea();
    }


    tmp<volScalarField> permeability() const
    {
        return permeabilityModelPtr_->permeability();
    }

    tmp<volScalarField> inversePermeability() const
    {
        return permeabilityModelPtr_->inversePermeability();
    }

    void updatePermeability()
    {
        return permeabilityModelPtr_->updatePermeability();
    }
*/

/*
    const PtrList<volScalarField>& mineralsDistribution() const
    {
        return Ys_;
    }

    const volScalarField & porosity() const
    {
        return eps_;
    }

*/

    void updatePorosity();

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

    void coucou()
    {
        Info << "COUCOU " << nl <<endl;
    }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
