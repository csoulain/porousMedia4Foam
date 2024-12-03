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
//#include "speciesTable.H"
#include "runTimeSelectionTables.H"
#include "unsaturatedPorousAssemblage.H"
#include "fluidProperties.H"
#include "incompressiblePhase.H"
//#include "autoPtr.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class twoPhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class basicUnsaturatedGeochemicalModel : public fluidProperties, public unsaturatedPorousAssemblage
{


  protected:

    // Protected data

    const fvMesh &mesh_;

    dictionary geochemicalModelDict_;

    dictionary fluidPropertiesDict_;

    volScalarField dMinvdRho_;
    volScalarField dGasvdRho_;


    word phiaName_;
    word phibName_;
    const surfaceScalarField &phia_;
    const surfaceScalarField &phib_;


    word phiName_;

    const surfaceScalarField &phi_;

    word SbName_;
 
    const volScalarField &Sb_;



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
              const dictionary& dict,
              const volScalarField &Sb,
              const incompressiblePhase &phasea,
              const incompressiblePhase &phaseb
          ),
          (mesh, dict, Sb, phasea, phaseb)
      );

    // Selectors

        //- Return a reference to the selected viscosity model
        static autoPtr<basicUnsaturatedGeochemicalModel> New
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const volScalarField &Sb,
            const incompressiblePhase &phasea,
            const incompressiblePhase &phaseb
        );


  // Constructors

    //- Construct from components

    basicUnsaturatedGeochemicalModel
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const volScalarField &Sb,
        const incompressiblePhase &phasea,
        const incompressiblePhase &phaseb
    );

    //- Destructor
    virtual ~basicUnsaturatedGeochemicalModel()
    {}


    // Member Functions

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
