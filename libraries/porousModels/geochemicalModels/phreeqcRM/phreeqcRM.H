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
    phreeqcRM module

SourceFiles
    phreeqc.C

\*---------------------------------------------------------------------------*/

#ifndef phreeqcRM_H
#define phreeqcRM_H

//#include "dictionary.H"
//#include "volFieldsFwd.H"
#include "volFields.H"
#include "speciesTable.H"
#include "runTimeSelectionTables.H"
#include "basicGeochemicalModel.H"
//#include "autoPtr.H"
#include "PhreeqcRM.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace geochemicalModels
{

/*---------------------------------------------------------------------------*\
                      Class twoPhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class phreeqcRM : public basicGeochemicalModel
{




  protected:

    // Protected data
    dictionary geochemicalModelDict_;

    dictionary  phreeqcDict_;

    Switch activateUpdatePorosity_ ;

    Switch setComponentH2O_ ;

    scalar numThreads_; // for openmp. Default 1.
    label nthread_;
    int nxyz_;

    scalar strangSteps_;

    word phreeqcInputFile_;
    word phreeqcDataBase_;

  //  speciesTable mineralList_;

    PtrList<dictionary> mineralSubDict_;

    PtrList<Switch> activatePhaseEquilibrium_ ;

    PtrList<dimensionedScalar> Vm_;

    volScalarField linkDomainZonesToPhreeqcSol_;

    Switch cvODE;
    double cvODETol_;
    Switch use_SNIA; 	//seq. non-iterative


    word UName_;
    const volVectorField & U_;


    volScalarField pH_;

    const word densitymodelType_;
    
    volScalarField rho_;

//    volScalarField Aee_;

    PhreeqcRM  phreeqc_;

    IRM_RESULT status;

  // Private Member Functions

    std::string generateEqPhasesInputString();

    std::string generateKineticsInputString();

    void updateKineticsParameters();

    void initializeMineralDistribution();

    void initializeFluidComposition();

    void activateSelectedOutput();

    void updatepH();

    void updatePorosityPhreeqc();

  // Private Member Functions

    //- Disallow copy construct
    phreeqcRM(const phreeqcRM&);


  public:

    //- Runtime type information
      TypeName("phreeqcRM");


    // Declare run-time constructor selection table

  // Constructors

    //- Construct from components

    phreeqcRM
    (
        const fvMesh& mesh,
        const dictionary& dict
    );

    //- Destructor
    virtual ~phreeqcRM()
    {}


    // Member Functions

    virtual void updateFluidComposition();

    virtual void updateMineralDistribution();

    virtual void updatePorosity();

    virtual tmp<volScalarField> rho() const;

    virtual void updateDensity();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace geochemicalModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
