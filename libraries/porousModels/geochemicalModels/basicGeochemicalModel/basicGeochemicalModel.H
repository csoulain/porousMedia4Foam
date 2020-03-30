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
    Foam::basicGeochemicalModel

Description
    basicGeochemicalModel module

SourceFiles
    basicGeochemicalModel.C

Authors
    17/02/2020 : C. Soulaine


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

    PtrList<volScalarField> dMs_;

    volScalarField inertMineral_;

    volScalarField eps_;

    PtrList<porousModel> porousMedia_;

//    tmp<volScalarField> Ae_;
//    volScalarField Ae_;

    autoPtr<absolutePermeabilityModel> absolutePermeabilityModelPtr_;

    autoPtr<dispersionModel> dispersionModelPtr_;


    word phiName_;

    const surfaceScalarField &phi_;

    // concentration fields
    PtrList<volScalarField> Y_;


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


    tmp<volScalarField> effectiveDispersion() const
    {
        return dispersionModelPtr_->effectiveDispersion();
    }

    void updateDispersion()
    {
        dispersionModelPtr_->updateDispersion();
    }



    const PtrList<volScalarField> dMs()
    {
        return dMs_;
    }

    virtual void updateFluidComposition() = 0;

    virtual void updateMineralDistribution() = 0;


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
